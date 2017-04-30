/*
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
  Adaptive controller by Ryan Beall

  See
  http://naira-hovakimyan.mechse.illinois.edu/l1-adaptive-control-tutorials/
  for mathematical basis
 */

#include <AP_HAL/AP_HAL.h>
#include <GCS_MAVLink/GCS.h>
#include "ADAP_Control.h"

extern const AP_HAL::HAL& hal;

const AP_Param::GroupInfo ADAP_Control::var_info[] = {
	// adaptive control parameters
	AP_GROUPINFO_FLAGS("CH", 1, ADAP_Control, enable_chan, 0, AP_PARAM_FLAG_ENABLE),
	AP_GROUPINFO("AL",       2, ADAP_Control, alpha, 24),
	AP_GROUPINFO("GAMT",     3, ADAP_Control, gamma_theta, 1000),
    AP_GROUPINFO("GAMW",     4, ADAP_Control, gamma_omega, 1000),
    AP_GROUPINFO("GAMS",     5, ADAP_Control, gamma_sigma, 1000),
	AP_GROUPINFO("THEU",     6, ADAP_Control, theta_max, 2.0),
	AP_GROUPINFO("THL",      7, ADAP_Control, theta_min, 0.5),
	AP_GROUPINFO("THE",      8, ADAP_Control, theta_epsilon, 5925),
    AP_GROUPINFO("OMU",      9, ADAP_Control, omega_max, 2.0),
	AP_GROUPINFO("OML",     10, ADAP_Control, omega_min, 0.5),
	AP_GROUPINFO("OME",     11, ADAP_Control, omega_epsilon, 5925),
    AP_GROUPINFO("SIGU",    12, ADAP_Control, sigma_max, 0.1),
	AP_GROUPINFO("SIGL",    13, ADAP_Control, sigma_min, -0.1),
	AP_GROUPINFO("SIGE",    14, ADAP_Control, sigma_epsilon, 203),
	AP_GROUPINFO("W0",      15, ADAP_Control, w0, 25),
    AP_GROUPINFO("K",       16, ADAP_Control, k, 0.1),
	AP_GROUPINFO("KG",      17, ADAP_Control, kg, 1.0),
    
	AP_GROUPEND
};

/*
  return true when enabled
 */
bool ADAP_Control::enabled(void) const
{
    return (enable_chan > 0 && hal.rcin->read(enable_chan-1) >= 1700);
}


/*
  reset to startup values
*/
void ADAP_Control::reset(uint16_t loop_rate_hz)
{
    x_m = x;       
    u = 0.0;
    u_lowpass = 0.0;
    u_sp = 0.0;
    theta = 1.0;
    omega = 1.0;
    sigma = 0.0;
    integrator = 0.0;

    float cutoff_hz = w0/(2*M_PI); //convert cutoff freq from rad/s to hz
    filter.set_cutoff_frequency(loop_rate_hz, cutoff_hz);
    filter.reset();
}

/*
  adaptive control update. Given a target rate in radians/second and a
  current sensor rate on the same axis in radians/second, return an
  actuator value from -1 to 1
 */
float ADAP_Control::update(uint16_t loop_rate_hz, float target_rate, float sensor_rate)
{
    float dt;

    x = sensor_rate;
    r = target_rate;

    //reset reference model at initialization
    uint64_t now = AP_HAL::micros64();
    if (last_run_us == 0 || now - last_run_us > 200000UL) {
        reset(loop_rate_hz);
        last_run_us = now;
        return 0;
    }

    dt = (now - last_run_us) * 1.0e-6;
    last_run_us = now;

    // u (controller output to plant)
    eta = theta*x + omega*u_lowpass + sigma;
    u_sp = eta;

    float out = constrain_float(eta-(kg*r),-radians(90)/dt, radians(90)/dt);

    //kD(s)
    integrator += dt*(out);  
    u = constrain_float(-k*integrator,-radians(45),radians(45)); // C(s)= wk/(s+wk) -> k sets the first order low pass response

    // Additional cascaded second order low pass filter (strictly propper)
    u_lowpass = filter.apply(u); 
    
    // State Predictor (first order single pole recursive filter)
    float alpha_filt = exp(-alpha*dt); //alpha in rad/s 
    alpha_filt = constrain_float(alpha_filt, 0.0, 1.0);
    float beta_filt = 1-alpha_filt;  

    x_m = alpha_filt*x_m + beta_filt*(u_sp);
       
    x_error = x_m - x;

    // Constrain error to +-300 deg/s
    x_error = constrain_float(x_error,-radians(300), radians(300));

    // Calculate solution to Lyapunov 1x1 matrix
    float Pb = 1/(2*alpha);

    // Projection Operator
    theta_dot = projection_operator(theta,-gamma_theta*x_error*Pb*x,theta_epsilon,theta_max,theta_min);
    omega_dot = projection_operator(omega,-gamma_omega*x_error*Pb*u_lowpass,omega_epsilon,omega_max,omega_min);
    sigma_dot = projection_operator(sigma,-gamma_sigma*x_error*Pb,sigma_epsilon,sigma_max,sigma_min);
			    
    // Parameter Update using Trapezoidal integration                 
    theta += dt*(theta_dot);
    omega += dt*(omega_dot);
    sigma += dt*(sigma_dot);

    theta = constrain_float(theta, theta_min, theta_max);
    omega = constrain_float(omega, omega_min, omega_max);
    sigma = constrain_float(sigma, sigma_min, sigma_max);
     
    // for ADAP_TUNING message
    theta_dot = theta_dot;
    omega_dot = omega_dot;
    sigma_dot = sigma_dot;
    x_error = x_error;

    DataFlash_Class::instance()->Log_Write(log_msg_name, "TimeUS,Dt,Atheta,Aomega,Asigma,Aeta,Axm,Ax,Ar,Axerr,AuL", "Qffffffffff",
                                           now,
                                           dt,
                                           theta, 
                                           omega,
                                           sigma,
                                           eta,
                                           degrees(x_m),
                                           degrees(x),
                                           degrees(r),
                                           degrees(x_error),
                                           degrees(u_lowpass));
    if (pid_info) {
        pid_info->P = theta;
        pid_info->I = sigma;
        pid_info->FF = x_m;
        pid_info->D = x_error;
        pid_info->desired = r;
    }
    
    return constrain_float(u_lowpass/radians(45), -1, 1);
}

float ADAP_Control::projection_operator(float Theta, float y, float epsilon, float th_max, float th_min) const
{
	
	// Calculate convex function
	// Nominal un-saturated value is above zero line on a parabolic curve
	// Steepness of curve is set by epsilon
	float f_diff2 = (th_max-th_min)*(th_max-th_min);
	float f_theta = (-4*(th_min - Theta) * (th_max - Theta))/(epsilon*f_diff2);
	float f_theta_dot = (4*(th_min + th_max - (2*Theta)))/(epsilon*f_diff2);	

	float projection_out = y;

	if (f_theta <= 0 && f_theta_dot*y < 0)
	{
		projection_out = y*(f_theta+1); //y-(y*(-1*f_theta));
	}	
 
   	return projection_out;
}

/*
  send ADAP_TUNING message
*/
void ADAP_Control::adaptive_tuning_send(mavlink_channel_t chan, uint8_t axis)
{
	if (!enabled() || !HAVE_PAYLOAD_SPACE(chan, ADAP_TUNING)) {
        return;
    }
    mavlink_msg_adap_tuning_send(chan, axis, 
                                 r,
                                 x,
                                 x_error,
                                 theta,
                                 omega,
                                 sigma,
                                 theta_dot,
                                 omega_dot,
                                 sigma_dot,
                                 x_m,	
                                 u_lowpass,
                                 u);
}

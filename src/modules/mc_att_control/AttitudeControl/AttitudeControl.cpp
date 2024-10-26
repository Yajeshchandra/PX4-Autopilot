/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file AttitudeControl.cpp
 */

#include <AttitudeControl.hpp>
#include <px4_platform_common/log.h>

#include <mathlib/math/Functions.hpp>
#include <matrix/Quaternion.hpp>


using namespace matrix;

void AttitudeControl::setProportionalGain(const matrix::Vector3f &proportional_gain, const float yaw_weight)
{
	_proportional_gain = proportional_gain;
	_yaw_w = math::constrain(yaw_weight, 0.f, 1.f);

	// compensate for the effect of the yaw weight rescaling the output
	if (_yaw_w > 1e-4f) {
		_proportional_gain(2) /= _yaw_w;
	}
}

void AttitudeControl::setPrimaryAxis(const Vector3f &axis, const float max_thrust)
{
	// Normalize the primary axis (n vector in paper)
	_primary_axis = axis.normalized();
	_max_thrust = max_thrust;

	// Equation 10: Thrust constraint check
	// -n_B_z >= m||g||/T_max
	if (-_primary_axis(2) * _max_thrust < 9.81f) {
		// Cannot maintain hover - default to level primary axis
		_primary_axis = Vector3f(0.0f, 0.0f, -1.0f);
	}
}


matrix::Vector3f AttitudeControl::update(const matrix::Quatf &q, const matrix::Vector3f &angular_velocity) const
{
	// Added angular_velocity as input since we need accurate yaw rate

	if (!_motor_failure) {
		return updateNormal(q);

	} else {
		// PX4_WARN("Motor failure detected, switching to degraded attitude control");
		return updateDegraded(q, angular_velocity);
	}
}

matrix::Vector3f AttitudeControl::updateNormal(const matrix::Quatf &q) const
{
	Quatf qd = _attitude_setpoint_q;

	// calculate reduced desired attitude neglecting vehicle's yaw to prioritize roll and pitch
	const Vector3f e_z = q.dcm_z();
	const Vector3f e_z_d = qd.dcm_z();
	Quatf qd_red(e_z, e_z_d);

	if (fabsf(qd_red(1)) > (1.f - 1e-5f) || fabsf(qd_red(2)) > (1.f - 1e-5f)) {
		// In the infinitesimal corner case where the vehicle and thrust have the completely opposite direction,
		// full attitude control anyways generates no yaw input and directly takes the combination of
		// roll and pitch leading to the correct desired yaw. Ignoring this case would still be totally safe and stable.
		qd_red = qd;

	} else {
		// transform rotation from current to desired thrust vector into a world frame reduced desired attitude
		qd_red *= q;
	}



	// mix full and reduced desired attitude
	matrix::Quatf q_mix = qd_red.inversed() * qd;
	q_mix.canonicalize();
	// catch numerical problems with the domain of acosf and asinf
	q_mix(0) = math::constrain(q_mix(0), -1.f, 1.f);
	q_mix(3) = math::constrain(q_mix(3), -1.f, 1.f);
	qd = qd_red * matrix::Quatf(cosf(_yaw_w * acosf(q_mix(0))), 0, 0, sinf(_yaw_w * asinf(q_mix(3))));

	// quaternion attitude control law, qe is rotation from q to qd
	const matrix::Quatf qe = q.inversed() * qd;

	// using sin(alpha/2) scaled rotation axis as attitude error (see quaternion definition by axis angle)
	// also taking care of the antipodal unit quaternion ambiguity
	const matrix::Vector3f eq = 2.f * qe.canonical().imag();

	// calculate angular rates setpoint
	matrix::Vector3f rate_setpoint = eq.emult(_proportional_gain);

	// Feed forward the yaw setpoint rate.
	// yawspeed_setpoint is the feed forward commanded rotation around the world z-axis,
	// but we need to apply it in the body frame (because _rates_sp is expressed in the body frame).
	// Therefore we infer the world z-axis (expressed in the body frame) by taking the last column of R.transposed (== q.inversed)
	// and multiply it by the yaw setpoint rate (yawspeed_setpoint).
	// This yields a vector representing the commanded rotatation around the world z-axis expressed in the body frame
	// such that it can be added to the rates setpoint.
	if (std::isfinite(_yawspeed_setpoint)) {
		rate_setpoint += q.inversed().dcm_z() * _yawspeed_setpoint;
	}

	// limit rates
	for (int i = 0; i < 3; i++) {
		rate_setpoint(i) = math::constrain(rate_setpoint(i), -_rate_limit(i), _rate_limit(i));
	}

	return rate_setpoint;
}

matrix::Vector3f AttitudeControl::updateDegraded(const Quatf &q, const Vector3f &angular_velocity) const
{
	// Implementation of Section III.B "Primary-Axis Attitude Control Loop"

	// Get current primary axis direction in inertial frame
	// n_I = L_IB * n_B (rotation of primary axis to inertial frame)
	matrix::Vector3f n_I = q.rotateVector(_primary_axis);
	(void)n_I;

	// Equation 11: Calculate n_des
	// n_des = (a_des - g)/||a_des - g||
	matrix::Vector3f g(0.0f, 0.0f, 9.81f);
	matrix::Vector3f a_des = _attitude_setpoint_q.dcm_z() * (_max_thrust / 9.81f);  // Scale by thrust
	matrix::Vector3f n_des = (a_des - g);
	n_des = n_des.normalized();

	// Equation 12: Express n_des in body frame
	// n_B_des = L_BI * n_I_des
	matrix::Vector3f n_des_B = q.inversed().rotateVector(n_des);

	// Equation 13: Extract components [h1 h2 h3]
	float h1 = n_des_B(0);
	float h2 = n_des_B(1);
	float h3 = n_des_B(2);

	// Get current yaw rate from angular velocity input
	// Paper uses r directly from gyro measurement
	float r = angular_velocity(2);  // Using z-axis angular velocity as yaw rate

	// Calculate n_des_dot in inertial frame for equation 15
	matrix::Vector3f n_des_dot = (n_des - _last_n_des) / _dt;
	matrix::Vector3f n_des_dot_B = q.inversed().rotateVector(n_des_dot);
	matrix::Vector2f n_hat_dot(n_des_dot_B(0), n_des_dot_B(1));  // Only need x,y components

	// Equation 17: Virtual input v_out
	// v_out = [n_dot_Bx + kx(n_Bx - h1); n_dot_By + ky(n_By - h2)]
	matrix::Vector2f v_out;
	v_out(0) = _k_primary(0) * (_primary_axis(0) - h1);  // n_dot_Bx = 0 (primary axis fixed)
	v_out(1) = _k_primary(1) * (_primary_axis(1) - h2);  // n_dot_By = 0

	// Equation 16: Calculate desired angular rates p_des and q_des
	// [p_des; q_des] = [0 1/h3; -1/h3 0] * (v_out - [h2; -h1]*r - n_hat_dot)
	matrix::Vector3f rate_setpoint;

	if (fabsf(h3) > 1e-4f) {  // Avoid division by small numbers
		// p_des = (-h3 * v_out_x - h2 * r - n_hat_dot_x) / h3
		rate_setpoint(0) = (-h3 * v_out(0) - h2 * r - n_hat_dot(0)) / h3;
		// q_des = (h3 * v_out_y + h1 * r - n_hat_dot_y) / h3
		rate_setpoint(1) = (h3 * v_out(1) + h1 * r - n_hat_dot(1)) / h3;

	} else {
		rate_setpoint(0) = 0.0f;
		rate_setpoint(1) = 0.0f;
	}

	// Keep measured yaw rate as discussed in Section III.B
	rate_setpoint(2) = r;

	// Apply rate limits
	for (int i = 0; i < 3; i++) {
		rate_setpoint(i) = math::constrain(rate_setpoint(i), -_rate_limit(i), _rate_limit(i));
	}

	_last_n_des = n_des;  // Store for next derivative calculation

	return rate_setpoint;
}

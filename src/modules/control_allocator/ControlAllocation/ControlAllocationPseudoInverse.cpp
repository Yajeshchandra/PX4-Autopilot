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
 * @file ControlAllocationPseudoInverse.hpp
 *
 * Simple Control Allocation Algorithm
 *
 * @author Julien Lecoeur <julien.lecoeur@gmail.com>
 */

#include "ControlAllocationPseudoInverse.hpp"


void
ControlAllocationPseudoInverse::setEffectivenessMatrix(
	const matrix::Matrix<float, ControlAllocation::NUM_AXES, ControlAllocation::NUM_ACTUATORS> &effectiveness,
	const ActuatorVector &actuator_trim, const ActuatorVector &linearization_point, int num_actuators,
	bool update_normalization_scale)
{
	ControlAllocation::setEffectivenessMatrix(effectiveness, actuator_trim, linearization_point, num_actuators,
			update_normalization_scale);
	_mix_update_needed = true;
	_normalization_needs_update = update_normalization_scale;
}


void
ControlAllocationPseudoInverse::setMotorFailure(int motor_idx, bool failed)
{
    if (motor_idx >= 0 && motor_idx < 4) {
        _motor_failure[motor_idx] = failed;
        _is_any_motor_failed = false;
        for (int i = 0; i < 4; i++) {
            if (_motor_failure[i]) {
                _is_any_motor_failed = true;
                break;
            }
        }
        _mix_update_needed = true;
        computePrimaryAxis();
    }
}

void
ControlAllocationPseudoInverse::computePrimaryAxis()
{
    if (_is_any_motor_failed) {
        // Set primary axis based on which motor failed
        // This follows the paper's approach of defining a primary rotation axis
        if (_motor_failure[3]) { // Left back motor (as in paper)
            _primary_axis = matrix::Vector3f(0.1f, 0.1f, -1.0f); // Initial values from paper
        } else if (_motor_failure[0]) {
            _primary_axis = matrix::Vector3f(-0.1f, -0.1f, -1.0f);
        } else if (_motor_failure[1]) {
            _primary_axis = matrix::Vector3f(-0.1f, 0.1f, -1.0f);
        } else if (_motor_failure[2]) {
            _primary_axis = matrix::Vector3f(0.1f, -0.1f, -1.0f);
        }
        _primary_axis.normalize();
    }
}

void
ControlAllocationPseudoInverse::updateEffectivenessForFailure()
{
    if (_is_any_motor_failed) {
        // Modify effectiveness matrix by zeroing out failed motor's column
        for (int i = 0; i < 4; i++) {
            if (_motor_failure[i]) {
                for (int j = 0; j < NUM_AXES; j++) {
                    _effectiveness(j, i) = 0.0f;
                }
            }
        }
    }
}

void
ControlAllocationPseudoInverse::allocate()
{
    // Update effectiveness matrix if motor failure state changed
    if (_mix_update_needed) {
        updateEffectivenessForFailure();
    }

    // Compute new gains if needed
    updatePseudoInverse();

    _prev_actuator_sp = _actuator_sp;

    if (!_is_any_motor_failed) {
        // Normal allocation for no failure case
        _actuator_sp = _actuator_trim + _mix * (_control_sp - _control_trim);
    } else {
        // Modified control allocation for failed state
        // Following paper's approach of primary axis control

        // Calculate desired thrust vector
        // Calculate desired thrust vector
        matrix::Vector3f thrust_desired(_control_sp(0), _control_sp(1), _control_sp(2));
        float thrust_magnitude = thrust_desired.length();

        // Project thrust onto primary axis
        float primary_axis_length_squared = _primary_axis.dot(_primary_axis);
        matrix::Vector3f thrust_primary = _primary_axis *
            (thrust_desired.dot(_primary_axis) / primary_axis_length_squared);

        // Calculate control inputs considering spinning motion
        matrix::Vector<float, NUM_AXES> modified_control = _control_sp;
        modified_control(0) = thrust_primary(0);
        modified_control(1) = thrust_primary(1);
        modified_control(2) = thrust_primary(2);


        // Allocate considering spinning motion
        _actuator_sp = _actuator_trim + _mix * (modified_control - _control_trim);

        // Scale outputs to maintain thrust while spinning
        float scale = thrust_magnitude / thrust_primary.length();
        for (int i = 0; i < _num_actuators; i++) {
            if (!_motor_failure[i]) {
                _actuator_sp(i) *= scale;
            } else {
                _actuator_sp(i) = 0.0f; // Ensure failed motor stays off
            }
        }
    }
}

void
ControlAllocationPseudoInverse::updatePseudoInverse()
{
    if (_mix_update_needed) {
        // Update effectiveness matrix for any motor failures
        updateEffectivenessForFailure();

        // Calculate pseudo-inverse with modified effectiveness matrix
        matrix::geninv(_effectiveness, _mix);

        if (_normalization_needs_update && !_had_actuator_failure) {
            updateControlAllocationMatrixScale();
            _normalization_needs_update = false;
        }

        normalizeControlAllocationMatrix();
        _mix_update_needed = false;
    }
}

void
ControlAllocationPseudoInverse::updateControlAllocationMatrixScale()
{
	// Same scale on roll and pitch
	if (_normalize_rpy) {

		int num_non_zero_roll_torque = 0;
		int num_non_zero_pitch_torque = 0;

		for (int i = 0; i < _num_actuators; i++) {

			if (fabsf(_mix(i, 0)) > 1e-3f) {
				++num_non_zero_roll_torque;
			}

			if (fabsf(_mix(i, 1)) > 1e-3f) {
				++num_non_zero_pitch_torque;
			}
		}

		float roll_norm_scale = 1.f;

		if (num_non_zero_roll_torque > 0) {
			roll_norm_scale = sqrtf(_mix.col(0).norm_squared() / (num_non_zero_roll_torque / 2.f));
		}

		float pitch_norm_scale = 1.f;

		if (num_non_zero_pitch_torque > 0) {
			pitch_norm_scale = sqrtf(_mix.col(1).norm_squared() / (num_non_zero_pitch_torque / 2.f));
		}

		_control_allocation_scale(0) = fmaxf(roll_norm_scale, pitch_norm_scale);
		_control_allocation_scale(1) = _control_allocation_scale(0);

		// Scale yaw separately
		_control_allocation_scale(2) = _mix.col(2).max();

	} else {
		_control_allocation_scale(0) = 1.f;
		_control_allocation_scale(1) = 1.f;
		_control_allocation_scale(2) = 1.f;
	}

	// Scale thrust by the sum of the individual thrust axes, and use the scaling for the Z axis if there's no actuators
	// (for tilted actuators)
	_control_allocation_scale(THRUST_Z) = 1.f;

	for (int axis_idx = 2; axis_idx >= 0; --axis_idx) {
		int num_non_zero_thrust = 0;
		float norm_sum = 0.f;

		for (int i = 0; i < _num_actuators; i++) {
			float norm = fabsf(_mix(i, 3 + axis_idx));
			norm_sum += norm;

			if (norm > FLT_EPSILON) {
				++num_non_zero_thrust;
			}
		}

		if (num_non_zero_thrust > 0) {
			_control_allocation_scale(3 + axis_idx) = norm_sum / num_non_zero_thrust;

		} else {
			_control_allocation_scale(3 + axis_idx) = _control_allocation_scale(THRUST_Z);
		}
	}
}

void
ControlAllocationPseudoInverse::normalizeControlAllocationMatrix()
{
	if (_control_allocation_scale(0) > FLT_EPSILON) {
		_mix.col(0) /= _control_allocation_scale(0);
		_mix.col(1) /= _control_allocation_scale(1);
	}

	if (_control_allocation_scale(2) > FLT_EPSILON) {
		_mix.col(2) /= _control_allocation_scale(2);
	}

	if (_control_allocation_scale(3) > FLT_EPSILON) {
		_mix.col(3) /= _control_allocation_scale(3);
		_mix.col(4) /= _control_allocation_scale(4);
		_mix.col(5) /= _control_allocation_scale(5);
	}

	// Set all the small elements to 0 to avoid issues
	// in the control allocation algorithms
	for (int i = 0; i < _num_actuators; i++) {
		for (int j = 0; j < NUM_AXES; j++) {
			if (fabsf(_mix(i, j)) < 1e-3f) {
				_mix(i, j) = 0.f;
			}
		}
	}
}

// void
// ControlAllocationPseudoInverse::allocate()
// {
// 	//Compute new gains if needed
// 	updatePseudoInverse();

// 	_prev_actuator_sp = _actuator_sp;

// 	// Allocate
// 	_actuator_sp = _actuator_trim + _mix * (_control_sp - _control_trim);
// }

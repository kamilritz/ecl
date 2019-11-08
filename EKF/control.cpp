/****************************************************************************
 *
 *   Copyright (c) 2015 Estimation and Control Library (ECL). All rights reserved.
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
 * 3. Neither the name ECL nor the names of its contributors may be
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
 * @file control.cpp
 * Control functions for ekf attitude and position estimator.
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include "../ecl.h"
#include "ekf.h"
#include <mathlib/mathlib.h>

void Ekf::controlFusionModes()
{
	// Store the status to enable change detection
	_control_status_prev.value = _control_status.value;

	// monitor the tilt alignment
	if (!_control_status.flags.tilt_align) {
		// whilst we are aligning the tilt, monitor the variances
		Vector3f angle_err_var_vec = calcRotVecVariances();

		// Once the tilt variances have reduced to equivalent of 3deg uncertainty, re-set the yaw and magnetic field states
		// and declare the tilt alignment complete
		if ((angle_err_var_vec(0) + angle_err_var_vec(1)) < sq(math::radians(3.0f))) {
			_control_status.flags.tilt_align = true;
			_control_status.flags.yaw_align = resetMagHeading(_mag_sample_delayed.mag);

			// send alignment status message to the console
			if (_control_status.flags.baro_hgt) {
				ECL_INFO("EKF aligned, (pressure height, IMU buf: %i, OBS buf: %i)", (int)_imu_buffer_length, (int)_obs_buffer_length);

			} else if (_control_status.flags.ev_hgt) {
				ECL_INFO("EKF aligned, (EV height, IMU buf: %i, OBS buf: %i)", (int)_imu_buffer_length, (int)_obs_buffer_length);

			} else if (_control_status.flags.gps_hgt) {
				ECL_INFO("EKF aligned, (GPS height, IMU buf: %i, OBS buf: %i)", (int)_imu_buffer_length, (int)_obs_buffer_length);

			} else if (_control_status.flags.rng_hgt) {
				ECL_INFO("EKF aligned, (range height, IMU buf: %i, OBS buf: %i)", (int)_imu_buffer_length, (int)_obs_buffer_length);
			} else {
				ECL_ERR("EKF aligned, (unknown height, IMU buf: %i, OBS buf: %i)", (int)_imu_buffer_length, (int)_obs_buffer_length);
			}

		}

	}

	// check for intermittent data (before pop_first_older_than)
	const baroSample &baro_init = _baro_buffer.get_newest();
	_baro_hgt_faulty = !((_time_last_imu - baro_init.time_us) < 2 * BARO_MAX_INTERVAL);

	const gpsSample &gps_init = _gps_buffer.get_newest();
	_gps_hgt_intermittent = !((_time_last_imu - gps_init.time_us) < 2 * GPS_MAX_INTERVAL);

	// check for arrival of new sensor data at the fusion time horizon
	_gps_data_ready = _gps_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_gps_sample_delayed);
	_mag_data_ready = _mag_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_mag_sample_delayed);

	if (_mag_data_ready) {
		// if enabled, use knowledge of theoretical magnetic field vector to calculate a synthetic magnetomter Z component value.
		// this is useful if there is a lot of interference on the sensor measurement.
		if (_params.synthesize_mag_z && (_params.mag_declination_source & MASK_USE_GEO_DECL) &&_NED_origin_initialised) {
			Vector3f mag_earth_pred = Dcmf(Eulerf(0, -_mag_inclination_gps, _mag_declination_gps)) * Vector3f(_mag_strength_gps, 0, 0);
			_mag_sample_delayed.mag(2) = calculate_synthetic_mag_z_measurement(_mag_sample_delayed.mag, mag_earth_pred);
			_control_status.flags.synthetic_mag_z = true;
		} else {
			_control_status.flags.synthetic_mag_z = false;
		}
	}

	_delta_time_baro_us = _baro_sample_delayed.time_us;
	_baro_data_ready = _baro_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_baro_sample_delayed);

	// if we have a new baro sample save the delta time between this sample and the last sample which is
	// used below for baro offset calculations
	if (_baro_data_ready) {
		_delta_time_baro_us = _baro_sample_delayed.time_us - _delta_time_baro_us;
	}

	// calculate 2,2 element of rotation matrix from sensor frame to earth frame
	// this is required for use of range finder and flow data
	_R_rng_to_earth_2_2 = _R_to_earth(2, 0) * _sin_tilt_rng + _R_to_earth(2, 2) * _cos_tilt_rng;

	// Get range data from buffer and check validity
	_range_data_ready = _range_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_range_sample_delayed);

	updateRangeDataValidity();

	if (_range_data_ready && _rng_hgt_valid) {
		// correct the range data for position offset relative to the IMU
		Vector3f pos_offset_body = _params.rng_pos_body - _params.imu_pos_body;
		Vector3f pos_offset_earth = _R_to_earth * pos_offset_body;
		_range_sample_delayed.rng += pos_offset_earth(2) / _R_rng_to_earth_2_2;
	}

	// We don't fuse flow data immediately because we have to wait for the mid integration point to fall behind the fusion time horizon.
	// This means we stop looking for new data until the old data has been fused.
	if (!_flow_data_ready) {
		_flow_data_ready = _flow_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_flow_sample_delayed)
				   && (_R_to_earth(2, 2) > _params.range_cos_max_tilt);
	}

	// check if we should fuse flow data for terrain estimation
	if (!_flow_for_terrain_data_ready && _flow_data_ready && _control_status.flags.in_air) {
		// only fuse flow for terrain if range data hasn't been fused for 5 seconds
		_flow_for_terrain_data_ready = (_time_last_imu - _time_last_hagl_fuse) > 5 * 1000 * 1000;
		// only fuse flow for terrain if the main filter is not fusing flow and we are using gps
		_flow_for_terrain_data_ready &= (!_control_status.flags.opt_flow && _control_status.flags.gps);
	}

	_ev_data_ready = _ext_vision_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_ev_sample_delayed);
	_tas_data_ready = _airspeed_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_airspeed_sample_delayed);
	_aux_vel_data_ready = _auxvel_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_auxvel_sample_delayed);

	// check for height sensor timeouts and reset and change sensor if necessary
	controlHeightSensorTimeouts();

	// control use of observations for aiding
	controlMagFusion();
	controlOpticalFlowFusion();
	controlGpsFusion();
	controlAirDataFusion();
	controlBetaFusion();
	controlDragFusion();
	controlHeightFusion();

	// Additional data odoemtery data from an external estimator can be fused.
	controlExternalVisionFusion();

	// Additional horizontal velocity data from an auxiliary sensor can be fused
	controlAuxVelFusion();

	// Fake position measurement for constraining drift when no other velocity or position measurements
	controlFakePosFusion();

	// check if we are no longer fusing measurements that directly constrain velocity drift
	update_deadreckoning_status();
}

void Ekf::controlExternalVisionFusion()
{
	// Check for new external vision data
	if (_ev_data_ready) {

		// if the ev data is not in a NED reference frame, then the transformation between EV and EKF navigation frames
		// needs to be calculated and the observations rotated into the EKF frame of reference
		if ((_params.fusion_mode & MASK_ROTATE_EV) && ((_params.fusion_mode & MASK_USE_EVPOS) || (_params.fusion_mode & MASK_USE_EVVEL)) && !_control_status.flags.ev_yaw) {
			// rotate EV measurements into the EKF Navigation frame
			calcExtVisRotMat();
		}

		// external vision aiding selection logic
		if (_control_status.flags.tilt_align && _control_status.flags.yaw_align) {

			// check for a external vision measurement that has fallen behind the fusion time horizon
			if ((_time_last_imu - _time_last_ext_vision) < (2 * EV_MAX_INTERVAL)) {
				// turn on use of external vision measurements for position
				if (_params.fusion_mode & MASK_USE_EVPOS && !_control_status.flags.ev_pos) {
					resetToGeneralHorizontalPosition();
					startEvPosFusion();
				}

				// turn on use of external vision measurements for velocity
				if (_params.fusion_mode & MASK_USE_EVVEL && !_control_status.flags.ev_vel) {
					resetToGeneralVelocity();
					startEvVelFusion();
				}

				if ((_params.fusion_mode & MASK_ROTATE_EV) && !(_params.fusion_mode & MASK_USE_EVYAW)
					&& !_ev_rot_mat_initialised)  {
					// Reset transformation between EV and EKF navigation frames when starting fusion
					resetExtVisRotMat();
					_ev_rot_mat_initialised = true;
					ECL_INFO_TIMESTAMPED("external vision aligned");
				}
			}
		}

		// external vision yaw aiding selection logic
		if (!_control_status.flags.gps && (_params.fusion_mode & MASK_USE_EVYAW) && !_control_status.flags.ev_yaw && _control_status.flags.tilt_align) {
			// don't start using EV data unless daa is arriving frequently
			if (_time_last_imu - _time_last_ext_vision < 2 * EV_MAX_INTERVAL) {
				// reset the yaw angle to the value from the observation quaternion
				// get the roll, pitch, yaw estimates from the quaternion states
				Quatf q_init(_state.quat_nominal);
				Eulerf euler_init(q_init);

				// get initial yaw from the observation quaternion
				const extVisionSample &ev_newest = _ext_vision_buffer.get_newest();
				Quatf q_obs(ev_newest.quat);
				Eulerf euler_obs(q_obs);
				euler_init(2) = euler_obs(2);

				// save a copy of the quaternion state for later use in calculating the amount of reset change
				Quatf quat_before_reset = _state.quat_nominal;

				// calculate initial quaternion states for the ekf
				_state.quat_nominal = Quatf(euler_init);
				uncorrelateQuatStates();

				// adjust the quaternion covariances estimated yaw error
				increaseQuatYawErrVariance(fmaxf(_ev_sample_delayed.angVar, 1.0e-4f));

				// calculate the amount that the quaternion has changed by
				_state_reset_status.quat_change = _state.quat_nominal * quat_before_reset.inversed();

				// add the reset amount to the output observer buffered data
				for (uint8_t i = 0; i < _output_buffer.get_length(); i++) {
					_output_buffer[i].quat_nominal = _state_reset_status.quat_change * _output_buffer[i].quat_nominal;
				}

				// apply the change in attitude quaternion to our newest quaternion estimate
				// which was already taken out from the output buffer
				_output_new.quat_nominal = _state_reset_status.quat_change * _output_new.quat_nominal;

				// capture the reset event
				_state_reset_status.quat_counter++;

				// flag the yaw as aligned
				_control_status.flags.yaw_align = true;

				startEvYawFusion();

			}
		}



		// determine if we should use the horizontal position observations
		if (_control_status.flags.ev_pos) {

			Vector3f ev_pos_obs_var{};
			Vector2f ev_pos_innov_gates{};

			// correct position and height for offset relative to IMU
			Vector3f pos_offset_body = _params.ev_pos_body - _params.imu_pos_body;
			Vector3f pos_offset_earth = _R_to_earth * pos_offset_body;
			_ev_sample_delayed.pos(0) -= pos_offset_earth(0);
			_ev_sample_delayed.pos(1) -= pos_offset_earth(1);
			_ev_sample_delayed.pos(2) -= pos_offset_earth(2);

			// Use an incremental position fusion method for EV position data if GPS is also used
			if (_params.fusion_mode & MASK_USE_GPS) {
				_fuse_hpos_as_odom = true;
			} else {
				_fuse_hpos_as_odom = false;
			}

			if (_fuse_hpos_as_odom) {
				if (!_hpos_prev_available) {
					// no previous observation available to calculate position change
					_fuse_pos = false;
					_hpos_prev_available = true;

				} else {
					// calculate the change in position since the last measurement
					Vector3f ev_delta_pos = _ev_sample_delayed.pos - _pos_meas_prev;

					// rotate measurement into body frame is required when fusing with GPS
					ev_delta_pos = _ev_rot_mat * ev_delta_pos;

					// use the change in position since the last measurement
					_ev_pos_innov(0) = _state.pos(0) - _hpos_pred_prev(0) - ev_delta_pos(0);
					_ev_pos_innov(1) = _state.pos(1) - _hpos_pred_prev(1) - ev_delta_pos(1);

					// observation 1-STD error, incremental pos observation is expected to have more uncertainty
					ev_pos_obs_var(0) = fmaxf(_ev_sample_delayed.posVar(0), 0.5f);
					ev_pos_obs_var(1) = fmaxf(_ev_sample_delayed.posVar(1), 0.5f);
				}

				// record observation and estimate for use next time
				_pos_meas_prev = _ev_sample_delayed.pos;
				_hpos_pred_prev(0) = _state.pos(0);
				_hpos_pred_prev(1) = _state.pos(1);

			} else {
				// use the absolute position
				Vector3f ev_pos_meas = _ev_sample_delayed.pos;
				if (_params.fusion_mode & MASK_ROTATE_EV) {
					ev_pos_meas = _ev_rot_mat * ev_pos_meas;
				}
				_ev_pos_innov(0) = _state.pos(0) - ev_pos_meas(0);
				_ev_pos_innov(1) = _state.pos(1) - ev_pos_meas(1);

				ev_pos_obs_var(0) = fmaxf(_ev_sample_delayed.posVar(0), 0.01f);
				ev_pos_obs_var(1) = fmaxf(_ev_sample_delayed.posVar(1), 0.01f);

				// check if we have been deadreckoning too long
				if (isTimedOut(_time_last_hor_pos_fuse, _params.reset_timeout_max)) {
					// reset velocity only if we have no another source of aiding constraining it
					if (isTimedOut(_time_last_of_fuse, (uint64_t)1E6) && isTimedOut(_time_last_hor_vel_fuse, (uint64_t)1E6)) {
						resetToGeneralVelocity();
					}
					resetToGeneralHorizontalPosition();
				}
			}

			// innovation gate size
			ev_pos_innov_gates(0) = fmaxf(_params.ev_pos_innov_gate, 1.0f);

			fuseHorizontalPosition(_ev_pos_innov, ev_pos_innov_gates, ev_pos_obs_var, _ev_pos_innov_var, _ev_pos_test_ratio);
		}

		// determine if we should use the velocity observations
		if (_control_status.flags.ev_vel) {

			Vector3f ev_vel_obs_var{};
			Vector2f ev_vel_innov_gates{};

			Vector3f vel_aligned{_ev_sample_delayed.vel};

			// rotate measurement into correct earth frame if required
			if (_params.fusion_mode & MASK_ROTATE_EV) {
				vel_aligned = _ev_rot_mat * _ev_sample_delayed.vel;
			}

			// correct velocity for offset relative to IMU
			Vector3f ang_rate = _imu_sample_delayed.delta_ang * (1.0f / _imu_sample_delayed.delta_ang_dt);
			Vector3f pos_offset_body = _params.ev_pos_body - _params.imu_pos_body;
			Vector3f vel_offset_body = cross_product(ang_rate, pos_offset_body);
			Vector3f vel_offset_earth = _R_to_earth * vel_offset_body;
			vel_aligned -= vel_offset_earth;

			_ev_vel_innov = _state.vel - vel_aligned;

			// check if we have been deadreckoning too long
			if (isTimedOut(_time_last_hor_pos_fuse, _params.reset_timeout_max)) {
				// reset velocity only if we have no another source of aiding constraining it
				if (isTimedOut(_time_last_of_fuse, (uint64_t)1E6) && isTimedOut(_time_last_hor_vel_fuse, (uint64_t)1E6)) {
					resetToGeneralVelocity();
				}
				resetToGeneralHorizontalPosition();
			}

			ev_vel_obs_var(0) = fmaxf(_ev_sample_delayed.velVar(0), 0.01f);
			ev_vel_obs_var(1) = fmaxf(_ev_sample_delayed.velVar(1), 0.01f);
			ev_vel_obs_var(2) = fmaxf(_ev_sample_delayed.velVar(2), 0.01f);

			ev_vel_innov_gates(0) = ev_vel_innov_gates(1) = fmaxf(_params.ev_vel_innov_gate, 1.0f);

			fuseHorizontalVelocity(_ev_vel_innov, ev_vel_innov_gates,ev_vel_obs_var, _ev_vel_innov_var, _ev_vel_test_ratio);
			fuseVerticalVelocity(_ev_vel_innov, ev_vel_innov_gates, ev_vel_obs_var, _ev_vel_innov_var, _ev_vel_test_ratio);
		}

		// determine if we should use the yaw observation
		if (_control_status.flags.ev_yaw) {
			fuseHeading();
		}

	} else if ((_control_status.flags.ev_pos || _control_status.flags.ev_vel)
		   && isTimedOut(_time_last_ext_vision,_params.reset_timeout_max)) {

		// Turn off EV fusion mode if no data has been received
		stopEvFusion();
		ECL_INFO_TIMESTAMPED("External Vision Data Stopped");

	}
}

void Ekf::controlOpticalFlowFusion()
{
	// Check if on ground motion is un-suitable for use of optical flow
	if (!_control_status.flags.in_air) {
		// When on ground check if the vehicle is being shaken or moved in a way that could cause a loss of navigation
		const float accel_norm = _accel_vec_filt.norm();

		const bool motion_is_excessive = ((accel_norm > (CONSTANTS_ONE_G * 1.5f)) // upper g limit
					    || (accel_norm < (CONSTANTS_ONE_G * 0.5f)) // lower g limit
					    || (_ang_rate_mag_filt > _flow_max_rate) // angular rate exceeds flow sensor limit
					    || (_R_to_earth(2,2) < cosf(math::radians(30.0f)))); // tilted excessively

		if (motion_is_excessive) {
			_time_bad_motion_us = _imu_sample_delayed.time_us;

		} else {
			_time_good_motion_us = _imu_sample_delayed.time_us;
		}

	} else {
		_time_bad_motion_us = 0;
		_time_good_motion_us = _imu_sample_delayed.time_us;
	}

	// Accumulate autopilot gyro data across the same time interval as the flow sensor
	_imu_del_ang_of += _imu_sample_delayed.delta_ang - _state.delta_ang_bias;
	_delta_time_of += _imu_sample_delayed.delta_ang_dt;

	// New optical flow data is available and is ready to be fused when the midpoint of the sample falls behind the fusion time horizon
	if (_flow_data_ready) {
		// Inhibit flow use if motion is un-suitable or we have good quality GPS
		// Apply hysteresis to prevent rapid mode switching
		float gps_err_norm_lim;
		if (_control_status.flags.opt_flow) {
			gps_err_norm_lim = 0.7f;
		} else {
			gps_err_norm_lim = 1.0f;
		}

		// Check if we are in-air and require optical flow to control position drift
		bool flow_required = _control_status.flags.in_air &&
				(_is_dead_reckoning // is doing inertial dead-reckoning so must constrain drift urgently
				  || (_control_status.flags.opt_flow && !_control_status.flags.gps && !_control_status.flags.ev_pos && !_control_status.flags.ev_vel) // is completely reliant on optical flow
				  || (_control_status.flags.gps && (_gps_error_norm > gps_err_norm_lim))); // is using GPS, but GPS is bad

		if (!_inhibit_flow_use && _control_status.flags.opt_flow) {
			// inhibit use of optical flow if motion is unsuitable and we are not reliant on it for flight navigation
			bool preflight_motion_not_ok = !_control_status.flags.in_air && ((_imu_sample_delayed.time_us - _time_good_motion_us) > (uint64_t)1E5);
			bool flight_motion_not_ok = _control_status.flags.in_air && !isRangeAidSuitable();
			if ((preflight_motion_not_ok || flight_motion_not_ok) && !flow_required) {
				_inhibit_flow_use = true;
			}
		} else if (_inhibit_flow_use && !_control_status.flags.opt_flow){
			// allow use of optical flow if motion is suitable or we are reliant on it for flight navigation
			bool preflight_motion_ok = !_control_status.flags.in_air && ((_imu_sample_delayed.time_us - _time_bad_motion_us) > (uint64_t)5E6);
			bool flight_motion_ok = _control_status.flags.in_air && isRangeAidSuitable();
			if (preflight_motion_ok || flight_motion_ok || flow_required) {
				_inhibit_flow_use = false;
			}
		}

		// Handle cases where we are using optical flow but are no longer able to because data is old
		// or its use has been inhibited.
		if (_control_status.flags.opt_flow) {
			if (_inhibit_flow_use) {
				stopFlowFusion();
				_time_last_of_fuse = 0;

			} else if ((_time_last_imu - _time_last_of_fuse) > (uint64_t)_params.reset_timeout_max) {
				stopFlowFusion();

			}
		}

		// optical flow fusion mode selection logic
		if ((_params.fusion_mode & MASK_USE_OF) // optical flow has been selected by the user
			&& !_control_status.flags.opt_flow // we are not yet using flow data
			&& _control_status.flags.tilt_align // we know our tilt attitude
			&& !_inhibit_flow_use
			&& isTerrainEstimateValid())
		{
			// If the heading is not aligned, reset the yaw and magnetic field states
			if (!_control_status.flags.yaw_align) {
				_control_status.flags.yaw_align = resetMagHeading(_mag_sample_delayed.mag);
			}

			// If the heading is valid and use is not inhibited , start using optical flow aiding
			if (_control_status.flags.yaw_align) {

			}

		}

		// handle the case when we have optical flow, are reliant on it, but have not been using it for an extended period
		if (_control_status.flags.opt_flow
		    && !_control_status.flags.gps
		    && !_control_status.flags.ev_pos
		    && !_control_status.flags.ev_vel) {

			bool do_reset = ((_time_last_imu - _time_last_of_fuse) > _params.reset_timeout_max);

			if (do_reset) {
				resetToGeneralVelocity();
				resetToGeneralHorizontalPosition();
			}
		}

		// Only fuse optical flow if valid body rate compensation data is available
		if (calcOptFlowBodyRateComp()) {

			bool flow_quality_good = (_flow_sample_delayed.quality >= _params.flow_qual_min);

			if (!flow_quality_good && !_control_status.flags.in_air) {
				// when on the ground with poor flow quality, assume zero ground relative velocity and LOS rate
				_flowRadXYcomp.zero();
			} else {
				// compensate for body motion to give a LOS rate
				_flowRadXYcomp(0) = _flow_sample_delayed.flowRadXY(0) - _flow_sample_delayed.gyroXYZ(0);
				_flowRadXYcomp(1) = _flow_sample_delayed.flowRadXY(1) - _flow_sample_delayed.gyroXYZ(1);
			}
		} else {
			// don't use this flow data and wait for the next data to arrive
			_flow_data_ready = false;
		}
	}

	// Wait until the midpoint of the flow sample has fallen behind the fusion time horizon
	if (_flow_data_ready && (_imu_sample_delayed.time_us > _flow_sample_delayed.time_us - uint32_t(1e6f * _flow_sample_delayed.dt) / 2)) {
		// Fuse optical flow LOS rate observations into the main filter only if height above ground has been updated recently
		// but use a relaxed time criteria to enable it to coast through bad range finder data
		if (_control_status.flags.opt_flow && ((_time_last_imu - _time_last_hagl_fuse) < (uint64_t)10e6)) {
			fuseOptFlow();
			_last_known_posNE(0) = _state.pos(0);
			_last_known_posNE(1) = _state.pos(1);
		}

		_flow_data_ready = false;
	}
}

void Ekf::controlGpsFusion()
{
	// Check for new GPS data that has fallen behind the fusion time horizon
	if (_gps_data_ready) {

		// GPS yaw aiding selection logic
		if ((_params.fusion_mode & MASK_USE_GPSYAW)
				&& ISFINITE(_gps_sample_delayed.yaw)
				&& _control_status.flags.tilt_align
				&& (!_control_status.flags.gps_yaw || !_control_status.flags.yaw_align)
				&& ((_time_last_imu - _time_last_gps) < (2 * GPS_MAX_INTERVAL))) {

			if (resetGpsAntYaw()) {
				// flag the yaw as aligned
				_control_status.flags.yaw_align = true;

				startGpsYawFusion();
			}
		}

		// fuse the yaw observation
		if (_control_status.flags.gps_yaw) {
			fuseGpsAntYaw();
		}

		// Determine if we should use GPS aiding for velocity and horizontal position
		// To start using GPS we need angular alignment completed, the local NED origin set and GPS data that has not failed checks recently
		bool gps_checks_passing = (_time_last_imu - _last_gps_fail_us > (uint64_t)5e6);
		bool gps_checks_failing = (_time_last_imu - _last_gps_pass_us > (uint64_t)5e6);

		if ((_params.fusion_mode & MASK_USE_GPS) && !_control_status.flags.gps) {
			if (_control_status.flags.tilt_align && _NED_origin_initialised && gps_checks_passing) {
				// If the heading is not aligned, reset the yaw and magnetic field states
				// Do not use external vision for yaw if using GPS because yaw needs to be
				// defined relative to an NED reference frame
				if (!_control_status.flags.yaw_align || _control_status.flags.ev_yaw || _mag_inhibit_yaw_reset_req) {
					_control_status.flags.yaw_align = resetMagHeading(_mag_sample_delayed.mag);
					if (_control_status.flags.yaw_align && _control_status.flags.ev_yaw) {
						stopEvYawFusion();
						ECL_INFO_TIMESTAMPED("stopped ev yaw fusion to start GPS fusion");
					}

					// Handle the special case where we have not been constraining yaw drift or learning yaw bias due
					// to assumed invalid mag field associated with indoor operation with a downwards looking flow sensor.
					if (_mag_inhibit_yaw_reset_req) {
						_mag_inhibit_yaw_reset_req = false;
						// Zero the yaw bias covariance and set the variance to the initial alignment uncertainty
						setDiag(P, 12, 12, sq(_params.switch_on_gyro_bias * FILTER_UPDATE_PERIOD_S));
					}
				}

				// If the heading is now aligned with North we can not start to fuse GPS position and velocities
				if (_control_status.flags.yaw_align) {
					startGpsFusion();
				}
			}
		}

		// Handle the case where we are using GPS and another source of aiding and GPS is failing checks
		if (_control_status.flags.gps  && gps_checks_failing && (_control_status.flags.opt_flow || _control_status.flags.ev_pos || _control_status.flags.ev_vel)) {
			stopGpsFusion();
			// Reset position state to external vision if we are going to use absolute values
			if (_control_status.flags.ev_pos && !(_params.fusion_mode & MASK_ROTATE_EV)) {
				resetToGeneralHorizontalPosition();
			}
			ECL_WARN_TIMESTAMPED("GPS data quality poor - stopping use");
		}

		// handle the case when we now have GPS, but have not been using it for an extended period
		if (_control_status.flags.gps) {
			// We are relying on aiding to constrain drift so after a specified time
			// with no aiding we need to do something
			bool do_reset = ((_time_last_imu - _time_last_hor_pos_fuse) > _params.reset_timeout_max)
					&& ((_time_last_imu - _time_last_delpos_fuse) > _params.reset_timeout_max)
					&& ((_time_last_imu - _time_last_hor_vel_fuse) > _params.reset_timeout_max)
					&& ((_time_last_imu - _time_last_of_fuse) > _params.reset_timeout_max);

			// We haven't had an absolute position fix for a longer time so need to do something
			do_reset = do_reset || ((_time_last_imu - _time_last_hor_pos_fuse) > (2 * _params.reset_timeout_max));

			if (do_reset) {
				// use GPS velocity data to check and correct yaw angle if a FW vehicle
				if (_control_status.flags.fixed_wing && _control_status.flags.in_air) {
					// if flying a fixed wing aircraft, do a complete reset that includes yaw
					_control_status.flags.mag_align_complete = realignYawGPS();
				}

				resetToGpsVelocity();
				resetToGpsHorizontalPosition();
				_velpos_reset_request = false;
				ECL_WARN_TIMESTAMPED("GPS fusion timeout - reset to GPS");

				// Reset the timeout counters
				_time_last_hor_pos_fuse = _time_last_imu;
				_time_last_hor_vel_fuse = _time_last_imu;

			}
		}

		// Only use GPS data for position and velocity aiding if enabled
		if (_control_status.flags.gps) {


			Vector2f gps_vel_innov_gates{}; // [horizontal vertical]
			Vector2f gps_pos_innov_gates{}; // [horizontal vertical]
			Vector3f gps_vel_obs_var{};
			Vector3f gps_pos_obs_var{};

			// correct velocity for offset relative to IMU
			Vector3f ang_rate = _imu_sample_delayed.delta_ang * (1.0f / _imu_sample_delayed.delta_ang_dt);
			Vector3f pos_offset_body = _params.gps_pos_body - _params.imu_pos_body;
			Vector3f vel_offset_body = cross_product(ang_rate, pos_offset_body);
			Vector3f vel_offset_earth = _R_to_earth * vel_offset_body;
			_gps_sample_delayed.vel -= vel_offset_earth;

			// correct position and height for offset relative to IMU
			Vector3f pos_offset_earth = _R_to_earth * pos_offset_body;
			_gps_sample_delayed.pos(0) -= pos_offset_earth(0);
			_gps_sample_delayed.pos(1) -= pos_offset_earth(1);
			_gps_sample_delayed.hgt += pos_offset_earth(2);

			// calculate observation process noise
			float lower_limit = fmaxf(_params.gps_pos_noise, 0.01f);

			if (_control_status.flags.opt_flow || _control_status.flags.ev_pos || _control_status.flags.ev_vel) {
				// if we are using other sources of aiding, then relax the upper observation
				// noise limit which prevents bad GPS perturbing the position estimate
				gps_pos_obs_var(0) = gps_pos_obs_var(1) = fmaxf(_gps_sample_delayed.hacc, lower_limit);

			} else {
				// if we are not using another source of aiding, then we are reliant on the GPS
				// observations to constrain attitude errors and must limit the observation noise value.
				float upper_limit = fmaxf(_params.pos_noaid_noise, lower_limit);
				gps_pos_obs_var(0) = gps_pos_obs_var(1) = math::constrain(_gps_sample_delayed.hacc, lower_limit, upper_limit);
			}

			gps_vel_obs_var(0) = gps_vel_obs_var(1) = gps_vel_obs_var(2) = sq(fmaxf(_gps_sample_delayed.sacc, _params.gps_vel_noise));

			// calculate innovations
			_gps_vel_innov(0) = _state.vel(0) - _gps_sample_delayed.vel(0);
			_gps_vel_innov(1) = _state.vel(1) - _gps_sample_delayed.vel(1);
			_gps_vel_innov(2) = _state.vel(2) - _gps_sample_delayed.vel(2);
			_gps_pos_innov(0) = _state.pos(0) - _gps_sample_delayed.pos(0);
			_gps_pos_innov(1) = _state.pos(1) - _gps_sample_delayed.pos(1);

			// set innovation gate size
			gps_pos_innov_gates(0) = fmaxf(_params.gps_pos_innov_gate, 1.0f);
			gps_pos_innov_gates(1) = fmaxf(_params.gps_vel_innov_gate, 1.0f);

			// fuse GPS measurement
			fuseHorizontalVelocity(_gps_vel_innov, gps_vel_innov_gates,gps_vel_obs_var, _gps_vel_innov_var, _gps_vel_test_ratio);
			fuseVerticalVelocity(_gps_vel_innov, gps_vel_innov_gates, gps_vel_obs_var, _gps_vel_innov_var, _gps_vel_test_ratio);
			fuseHorizontalPosition(_gps_pos_innov, gps_pos_innov_gates, gps_pos_obs_var, _gps_pos_innov_var, _gps_pos_test_ratio);
		}

	} else if (_control_status.flags.gps && isTimedOut(_gps_sample_delayed.time_us, (uint64_t)10e6)) {
		stopGpsFusion();
		ECL_WARN_TIMESTAMPED("GPS data stopped");
	}  else if (_control_status.flags.gps && isTimedOut(_gps_sample_delayed.time_us, (uint64_t)1e6) && (_control_status.flags.opt_flow || _control_status.flags.ev_pos || _control_status.flags.ev_vel)) {
		// Handle the case where we are fusing another position source along GPS,
		// stop waiting for GPS after 1 s of lost signal
		stopGpsFusion();
		ECL_WARN_TIMESTAMPED("GPS data stopped, using only EV or OF");
	}
}

void Ekf::controlHeightSensorTimeouts()
{
	/*
	 * Handle the case where we have not fused height measurements recently and
	 * uncertainty exceeds the max allowable. Reset using the best available height
	 * measurement source, continue using it after the reset and declare the current
	 * source failed if we have switched.
	*/

	// Check for IMU accelerometer vibration induced clipping as evidenced by the vertical innovations being positive and not stale.
	// Clipping causes the average accel reading to move towards zero which makes the INS think it is falling and produces positive vertical innovations
	float var_product_lim = sq(_params.vert_innov_test_lim) * sq(_params.vert_innov_test_lim);
	bool bad_vert_accel = (_control_status.flags.baro_hgt && // we can only run this check if vertical position and velocity observations are independent
			(sq(_gps_pos_innov(2) * fmaxf(fabsf(_gps_vel_innov(2)),fabsf(_ev_vel_innov(2)))) > var_product_lim * (_gps_pos_innov_var(2) * fmaxf(fabsf(_gps_vel_innov_var(2)),fabsf(_ev_vel_innov_var(2))))) && // vertical position and velocity sensors are in agreement that we have a significant error
			(_gps_vel_innov(2) > 0.0f || _ev_vel_innov(2) > 0.0f) && // positive innovation indicates that the inertial nav thinks it is falling
			((_imu_sample_delayed.time_us - _baro_sample_delayed.time_us) < 2 * BARO_MAX_INTERVAL) && // vertical position data is fresh
			((_imu_sample_delayed.time_us - _gps_sample_delayed.time_us) < 2 * GPS_MAX_INTERVAL)); // vertical velocity data is fresh

	// record time of last bad vert accel
	if (bad_vert_accel) {
		_time_bad_vert_accel =  _time_last_imu;

	} else {
		_time_good_vert_accel = _time_last_imu;
	}

	// declare a bad vertical acceleration measurement and make the declaration persist
	// for a minimum of 10 seconds
	if (_bad_vert_accel_detected) {
		_bad_vert_accel_detected = (_time_last_imu - _time_bad_vert_accel < BADACC_PROBATION);

	} else {
		_bad_vert_accel_detected = bad_vert_accel;
	}

	// check if height is continuously failing because of accel errors
	bool continuous_bad_accel_hgt = ((_time_last_imu - _time_good_vert_accel) > (unsigned)_params.bad_acc_reset_delay_us);

	// check if height has been inertial deadreckoning for too long
	bool hgt_fusion_timeout = isTimedOut(_time_last_hgt_fuse, (uint64_t)5e6);

	// reset the vertical position and velocity states
	if (hgt_fusion_timeout || continuous_bad_accel_hgt) {
		// boolean that indicates we will do a height reset
		bool reset_height = false;

		// handle the case where we are using baro for height
		if (_control_status.flags.baro_hgt) {
			// check if GPS height is available
			const gpsSample &gps_init = _gps_buffer.get_newest();
			bool gps_hgt_accurate = (gps_init.vacc < _params.req_vacc);

			const baroSample &baro_init = _baro_buffer.get_newest();
			bool baro_hgt_available = ((_time_last_imu - baro_init.time_us) < 2 * BARO_MAX_INTERVAL);

			// check for inertial sensing errors in the last 10 seconds
			bool prev_bad_vert_accel = (_time_last_imu - _time_bad_vert_accel < BADACC_PROBATION);

			// reset to GPS if adequate GPS data is available and the timeout cannot be blamed on IMU data
			bool reset_to_gps = !_gps_hgt_intermittent && gps_hgt_accurate && !prev_bad_vert_accel;

			// reset to GPS if GPS data is available and there is no Baro data
			reset_to_gps = reset_to_gps || (!_gps_hgt_intermittent && !baro_hgt_available);

			// reset to Baro if we are not doing a GPS reset and baro data is available
			bool reset_to_baro = !reset_to_gps && baro_hgt_available;

			if (reset_to_gps) {
				// set height sensor health
				_baro_hgt_faulty = true;

				// reset the height mode
				setControlGPSHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("baro hgt timeout - reset to GPS");

			} else if (reset_to_baro) {
				// set height sensor health
				_baro_hgt_faulty = false;

				// reset the height mode
				setControlBaroHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("baro hgt timeout - reset to baro");

			} else {
				// we have nothing we can reset to
				// deny a reset
				reset_height = false;

			}
		}

		// handle the case we are using GPS for height
		if (_control_status.flags.gps_hgt) {
			// check if GPS height is available
			const gpsSample &gps_init = _gps_buffer.get_newest();
			bool gps_hgt_accurate = (gps_init.vacc < _params.req_vacc);

			// check the baro height source for consistency and freshness
			const baroSample &baro_init = _baro_buffer.get_newest();
			bool baro_data_fresh = ((_time_last_imu - baro_init.time_us) < 2 * BARO_MAX_INTERVAL);
			float baro_innov = _state.pos(2) - (_hgt_sensor_offset - baro_init.hgt + _baro_hgt_offset);
			bool baro_data_consistent = fabsf(baro_innov) < (sq(_params.baro_noise) + P[9][9]) * sq(_params.baro_innov_gate);

			// if baro data is acceptable and GPS data is inaccurate, reset height to baro
			bool reset_to_baro = baro_data_consistent && baro_data_fresh && !_baro_hgt_faulty && !gps_hgt_accurate;

			// if GPS height is unavailable and baro data is available, reset height to baro
			reset_to_baro = reset_to_baro || (_gps_hgt_intermittent && baro_data_fresh);

			// if we cannot switch to baro and GPS data is available, reset height to GPS
			bool reset_to_gps = !reset_to_baro && !_gps_hgt_intermittent;

			if (reset_to_baro) {
				// set height sensor health
				_baro_hgt_faulty = false;

				// reset the height mode
				setControlBaroHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("gps hgt timeout - reset to baro");

			} else if (reset_to_gps) {
				// reset the height mode
				setControlGPSHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("gps hgt timeout - reset to GPS");

			} else {
				// we have nothing to reset to
				reset_height = false;

			}
		}

		// handle the case we are using range finder for height
		if (_control_status.flags.rng_hgt) {

			// check if baro data is available
			const baroSample &baro_init = _baro_buffer.get_newest();
			bool baro_data_available = ((_time_last_imu - baro_init.time_us) < 2 * BARO_MAX_INTERVAL);

			// reset to baro if we have no range data and baro data is available
			bool reset_to_baro = !_rng_hgt_valid && baro_data_available;

			if (_rng_hgt_valid) {

				// reset the height mode
				setControlRangeHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("rng hgt timeout - reset to rng hgt");

			} else if (reset_to_baro) {
				// set height sensor health
				_baro_hgt_faulty = false;

				// reset the height mode
				setControlBaroHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("rng hgt timeout - reset to baro");

			} else {
				// we have nothing to reset to
				reset_height = false;

			}
		}

		// handle the case where we are using external vision data for height
		if (_control_status.flags.ev_hgt) {
			// check if vision data is available
			const extVisionSample &ev_init = _ext_vision_buffer.get_newest();
			bool ev_data_available = ((_time_last_imu - ev_init.time_us) < 2 * EV_MAX_INTERVAL);

			// check if baro data is available
			const baroSample &baro_init = _baro_buffer.get_newest();
			bool baro_data_available = ((_time_last_imu - baro_init.time_us) < 2 * BARO_MAX_INTERVAL);

			// reset to baro if we have no vision data and baro data is available
			bool reset_to_baro = !ev_data_available && baro_data_available;

			// reset to ev data if it is available
			bool reset_to_ev = ev_data_available;

			if (reset_to_baro) {
				// set height sensor health
				_baro_hgt_faulty = false;

				// reset the height mode
				setControlBaroHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("ev hgt timeout - reset to baro");

			} else if (reset_to_ev) {
				// reset the height mode
				setControlEVHeight();

				// request a reset
				reset_height = true;
				ECL_WARN_TIMESTAMPED("ev hgt timeout - reset to ev hgt");

			} else {
				// we have nothing to reset to
				reset_height = false;

			}
		}

		// Reset vertical position and velocity states to the last measurement
		if (reset_height) {
			resetHeight();
			// Reset the timout timer
			_time_last_hgt_fuse = _time_last_imu;

		}

	}
}

void Ekf::controlHeightFusion()
{

	checkRangeAidSuitability();
	_range_aid_mode_selected = (_params.range_aid == 1) && isRangeAidSuitable();

	if (_params.vdist_sensor_type == VDIST_SENSOR_BARO) {

		if (_range_aid_mode_selected && _range_data_ready && _rng_hgt_valid) {
			setControlRangeHeight();
			_fuse_height = true;

			// we have just switched to using range finder, calculate height sensor offset such that current
			// measurement matches our current height estimate
			if (_control_status_prev.flags.rng_hgt != _control_status.flags.rng_hgt) {
				if (isTerrainEstimateValid()) {
					_hgt_sensor_offset = _terrain_vpos;

				} else {
					_hgt_sensor_offset = _R_rng_to_earth_2_2 * _range_sample_delayed.rng + _state.pos(2);
				}
			}

		} else if (!_range_aid_mode_selected && _baro_data_ready && !_baro_hgt_faulty) {
			setControlBaroHeight();
			_fuse_height = true;

			// we have just switched to using baro height, we don't need to set a height sensor offset
			// since we track a separate _baro_hgt_offset
			if (_control_status_prev.flags.baro_hgt != _control_status.flags.baro_hgt) {
				_hgt_sensor_offset = 0.0f;
			}

			// Turn off ground effect compensation if it times out
			if (_control_status.flags.gnd_effect) {
				if ((_time_last_imu - _time_last_gnd_effect_on > GNDEFFECT_TIMEOUT)) {

					_control_status.flags.gnd_effect = false;
				}
			}

		} else if (_control_status.flags.gps_hgt && _gps_data_ready && !_gps_hgt_intermittent) {
			// switch to gps if there was a reset to gps
			_fuse_height = true;

			// we have just switched to using gps height, calculate height sensor offset such that current
			// measurement matches our current height estimate
			if (_control_status_prev.flags.gps_hgt != _control_status.flags.gps_hgt) {
				_hgt_sensor_offset = _gps_sample_delayed.hgt - _gps_alt_ref + _state.pos(2);
			}
		}
	}

	// set the height data source to range if requested
	if ((_params.vdist_sensor_type == VDIST_SENSOR_RANGE) && _rng_hgt_valid) {
		setControlRangeHeight();
		_fuse_height = _range_data_ready;

		// we have just switched to using range finder, calculate height sensor offset such that current
		// measurement matches our current height estimate
		if (_control_status_prev.flags.rng_hgt != _control_status.flags.rng_hgt) {
			// use the parameter rng_gnd_clearance if on ground to avoid a noisy offset initialization (e.g. sonar)
			if (_control_status.flags.in_air && isTerrainEstimateValid()) {

				_hgt_sensor_offset = _terrain_vpos;

			} else if (_control_status.flags.in_air) {

				_hgt_sensor_offset = _R_rng_to_earth_2_2 * _range_sample_delayed.rng + _state.pos(2);

			} else {

				_hgt_sensor_offset = _params.rng_gnd_clearance;
			}
		}

	} else if ((_params.vdist_sensor_type == VDIST_SENSOR_RANGE) && _baro_data_ready && !_baro_hgt_faulty) {
		setControlBaroHeight();
		_fuse_height = true;

		// we have just switched to using baro height, we don't need to set a height sensor offset
		// since we track a separate _baro_hgt_offset
		if (_control_status_prev.flags.baro_hgt != _control_status.flags.baro_hgt) {
			_hgt_sensor_offset = 0.0f;
		}
	}

	// Determine if GPS should be used as the height source
	if (_params.vdist_sensor_type == VDIST_SENSOR_GPS) {

		if (_range_aid_mode_selected && _range_data_ready && _rng_hgt_valid) {
			setControlRangeHeight();
			_fuse_height = true;

			// we have just switched to using range finder, calculate height sensor offset such that current
			// measurement matches our current height estimate
			if (_control_status_prev.flags.rng_hgt != _control_status.flags.rng_hgt) {
				if (isTerrainEstimateValid()) {
					_hgt_sensor_offset = _terrain_vpos;

				} else {
					_hgt_sensor_offset = _R_rng_to_earth_2_2 * _range_sample_delayed.rng + _state.pos(2);
				}
			}

		} else if (!_range_aid_mode_selected && _gps_data_ready && !_gps_hgt_intermittent && _gps_checks_passed) {
			setControlGPSHeight();
			_fuse_height = true;

			// we have just switched to using gps height, calculate height sensor offset such that current
			// measurement matches our current height estimate
			if (_control_status_prev.flags.gps_hgt != _control_status.flags.gps_hgt) {
				_hgt_sensor_offset = _gps_sample_delayed.hgt - _gps_alt_ref + _state.pos(2);
			}

		} else if (_control_status.flags.baro_hgt && _baro_data_ready && !_baro_hgt_faulty) {
			// switch to baro if there was a reset to baro
			_fuse_height = true;

			// we have just switched to using baro height, we don't need to set a height sensor offset
			// since we track a separate _baro_hgt_offset
			if (_control_status_prev.flags.baro_hgt != _control_status.flags.baro_hgt) {
				_hgt_sensor_offset = 0.0f;
			}
		}
	}

	// Determine if we rely on EV height but switched to baro
	if (_params.vdist_sensor_type == VDIST_SENSOR_EV) {

		// don't start using EV data unless data is arriving frequently
		if (!_control_status.flags.ev_hgt && ((_time_last_imu - _time_last_ext_vision) < (2 * EV_MAX_INTERVAL))) {
			_fuse_height = true;
			setControlEVHeight();
			resetHeight();
		}

		if (_control_status.flags.baro_hgt && _baro_data_ready && !_baro_hgt_faulty) {
			// switch to baro if there was a reset to baro
			_fuse_height = true;

			// we have just switched to using baro height, we don't need to set a height sensor offset
			// since we track a separate _baro_hgt_offset
			if (_control_status_prev.flags.baro_hgt != _control_status.flags.baro_hgt) {
				_hgt_sensor_offset = 0.0f;
			}
		}
		// TODO: Add EV normal case here
		// determine if we should use the vertical position observation
		if (_control_status.flags.ev_hgt) {
			_fuse_height = true;
		}
	}

	// calculate a filtered offset between the baro origin and local NED origin if we are not using the baro as a height reference
	if (!_control_status.flags.baro_hgt && _baro_data_ready) {
		float local_time_step = 1e-6f * _delta_time_baro_us;
		local_time_step = math::constrain(local_time_step, 0.0f, 1.0f);

		// apply a 10 second first order low pass filter to baro offset
		float offset_rate_correction =  0.1f * (_baro_sample_delayed.hgt + _state.pos(
				2) - _baro_hgt_offset);
		_baro_hgt_offset += local_time_step * math::constrain(offset_rate_correction, -0.1f, 0.1f);
	}

	if ((_time_last_imu - _time_last_hgt_fuse) > 2 * RNG_MAX_INTERVAL && _control_status.flags.rng_hgt
	    && (!_range_data_ready || !_rng_hgt_valid)) {

		// If we are supposed to be using range finder data as the primary height sensor, have missed or rejected measurements
		// and are on the ground, then synthesise a measurement at the expected on ground value
		if (!_control_status.flags.in_air) {
			_range_sample_delayed.rng = _params.rng_gnd_clearance;
			_range_sample_delayed.time_us = _imu_sample_delayed.time_us;

		}

		_fuse_height = true;
	}

	if (_fuse_height) {

		Vector2f height_innov_gate{};
		Vector2f height_test_ratio{};
		Vector3f height_obs_var{};
		Vector3f height_innov_var{};
		Vector3f height_innov{};

		if (_control_status.flags.baro_hgt) {
			// vertical position innovation - baro measurement has opposite sign to earth z axis
			height_innov(2) = _state.pos(2) + _baro_sample_delayed.hgt - _baro_hgt_offset - _hgt_sensor_offset;
			// observation variance - user parameter defined
			height_obs_var(2) = sq(fmaxf(_params.baro_noise, 0.01f));

			// innovation gate size
			height_innov_gate(1) = fmaxf(_params.baro_innov_gate, 1.0f);

			// Compensate for positive static pressure transients (negative vertical position innovations)
			// caused by rotor wash ground interaction by applying a temporary deadzone to baro innovations.
			float deadzone_start = 0.0f;
			float deadzone_end = deadzone_start + _params.gnd_effect_deadzone;

			if (_control_status.flags.gnd_effect) {
				if (height_innov(2) < -deadzone_start) {
					if (height_innov(2) <= -deadzone_end) {
						height_innov(2) += deadzone_end;

					} else {
						height_innov(2) = -deadzone_start;
					}
				}
			}

		} else if (_control_status.flags.gps_hgt) {
			// vertical position innovation - gps measurement has opposite sign to earth z axis
			height_innov(2) = _state.pos(2) + _gps_sample_delayed.hgt - _gps_alt_ref - _hgt_sensor_offset;
			// observation variance - receiver defined and parameter limited
			// use scaled horizontal position accuracy assuming typical ratio of VDOP/HDOP
			float lower_limit = fmaxf(_params.gps_pos_noise, 0.01f);
			float upper_limit = fmaxf(_params.pos_noaid_noise, lower_limit);
			height_obs_var(2) = sq(1.5f * math::constrain(_gps_sample_delayed.vacc, lower_limit, upper_limit));
			// innovation gate size
			height_innov_gate(1) = fmaxf(_params.baro_innov_gate, 1.0f);

		} else if (_control_status.flags.rng_hgt && (_R_rng_to_earth_2_2 > _params.range_cos_max_tilt)) {
			// use range finder with tilt correction
			height_innov(2) = _state.pos(2) - (-math::max(_range_sample_delayed.rng * _R_rng_to_earth_2_2,
							 _params.rng_gnd_clearance)) - _hgt_sensor_offset;
			// observation variance - user parameter defined
			height_obs_var(2) = fmaxf((sq(_params.range_noise) + sq(_params.range_noise_scaler * _range_sample_delayed.rng)) * sq(_R_rng_to_earth_2_2), 0.01f);
			// innovation gate size
			height_innov_gate(1) = fmaxf(_params.range_innov_gate, 1.0f);

		} else if (_control_status.flags.ev_hgt) {
			// calculate the innovation assuming the external vision observation is in local NED frame
			height_innov(2) = _state.pos(2) - _ev_sample_delayed.pos(2);
			// observation variance - defined externally
			height_obs_var(2) = fmaxf(_ev_sample_delayed.posVar(2), 1e-4f);
			// innovation gate size
			height_innov_gate(1) = fmaxf(_params.ev_pos_innov_gate, 1.0f);
		}
		// fuse height inforamtion
		fuseVerticalPosition(height_innov,height_innov_gate,
				height_obs_var, height_innov_var,height_test_ratio);

		// This is a temporary hack until we do proper height sensor fusion
		_gps_pos_innov(2) = height_innov(2);
		_gps_pos_innov_var(2) = height_innov_var(2);
		_gps_pos_test_ratio(1) = height_test_ratio(1);
	}

}

void Ekf::checkRangeAidSuitability()
{
	const bool horz_vel_valid = _control_status.flags.gps
				    || _control_status.flags.ev_pos
				    || _control_status.flags.ev_vel
				    || _control_status.flags.opt_flow;

	if (_control_status.flags.in_air
	    && _rng_hgt_valid
	    && isTerrainEstimateValid()
	    && horz_vel_valid) {
		// check if we can use range finder measurements to estimate height, use hysteresis to avoid rapid switching
		// Note that the 0.7 coefficients and the innovation check are arbitrary values but work well in practice
		const bool is_in_range = _is_range_aid_suitable
					 ? (_terrain_vpos - _state.pos(2) < _params.max_hagl_for_range_aid)
					 : (_terrain_vpos - _state.pos(2) < _params.max_hagl_for_range_aid * 0.7f);

		const float ground_vel = sqrtf(_state.vel(0) * _state.vel(0) + _state.vel(1) * _state.vel(1));
		const bool is_below_max_speed = _is_range_aid_suitable
						? ground_vel < _params.max_vel_for_range_aid
						: ground_vel < _params.max_vel_for_range_aid * 0.7f;

		const bool is_hagl_stable = _is_range_aid_suitable
					    ? ((_hagl_innov * _hagl_innov / (sq(_params.range_aid_innov_gate) * _hagl_innov_var)) < 1.0f)
					    : ((_hagl_innov * _hagl_innov / (sq(_params.range_aid_innov_gate) * _hagl_innov_var)) < 0.01f);

		_is_range_aid_suitable = is_in_range && is_below_max_speed && is_hagl_stable;

	} else {
		_is_range_aid_suitable = false;
	}
}

void Ekf::controlAirDataFusion()
{
	// control activation and initialisation/reset of wind states required for airspeed fusion

	// If both airspeed and sideslip fusion have timed out and we are not using a drag observation model then we no longer have valid wind estimates
	bool airspeed_timed_out = ((_time_last_imu - _time_last_arsp_fuse) > (uint64_t)10e6);
	bool sideslip_timed_out = ((_time_last_imu - _time_last_beta_fuse) > (uint64_t)10e6);

	if (_control_status.flags.wind && airspeed_timed_out && sideslip_timed_out && !(_params.fusion_mode & MASK_USE_DRAG)) {
		_control_status.flags.wind = false;

	}

	if (_control_status.flags.fuse_aspd && airspeed_timed_out) {
		_control_status.flags.fuse_aspd = false;

	}

	// Always try to fuse airspeed data if available and we are in flight
	if (_tas_data_ready && _control_status.flags.in_air) {
		// always fuse airsped data if we are flying and data is present
		if (!_control_status.flags.fuse_aspd) {
			_control_status.flags.fuse_aspd = true;
		}

		// If starting wind state estimation, reset the wind states and covariances before fusing any data
		if (!_control_status.flags.wind) {
			// activate the wind states
			_control_status.flags.wind = true;
			// reset the timout timer to prevent repeated resets
			_time_last_arsp_fuse = _time_last_imu;
			_time_last_beta_fuse = _time_last_imu;
			// reset the wind speed states and corresponding covariances
			resetWindStates();
			resetWindCovariance();

		}

		fuseAirspeed();

	}
}

void Ekf::controlBetaFusion()
{
	// control activation and initialisation/reset of wind states required for synthetic sideslip fusion fusion

	// If both airspeed and sideslip fusion have timed out and we are not using a drag observation model then we no longer have valid wind estimates
	bool sideslip_timed_out = ((_time_last_imu - _time_last_beta_fuse) > (uint64_t)10e6);
	bool airspeed_timed_out = ((_time_last_imu - _time_last_arsp_fuse) > (uint64_t)10e6);

	if (_control_status.flags.wind && airspeed_timed_out && sideslip_timed_out && !(_params.fusion_mode & MASK_USE_DRAG)) {
		_control_status.flags.wind = false;
	}

	// Perform synthetic sideslip fusion when in-air and sideslip fuson had been enabled externally in addition to the following criteria:

	// Sufficient time has lapsed sice the last fusion
	bool beta_fusion_time_triggered = ((_time_last_imu - _time_last_beta_fuse) > _params.beta_avg_ft_us);

	if (beta_fusion_time_triggered && _control_status.flags.fuse_beta && _control_status.flags.in_air) {
		// If starting wind state estimation, reset the wind states and covariances before fusing any data
		if (!_control_status.flags.wind) {
			// activate the wind states
			_control_status.flags.wind = true;
			// reset the timeout timers to prevent repeated resets
			_time_last_beta_fuse = _time_last_imu;
			_time_last_arsp_fuse = _time_last_imu;
			// reset the wind speed states and corresponding covariances
			resetWindStates();
			resetWindCovariance();
		}

		fuseSideslip();
	}
}

void Ekf::controlDragFusion()
{
	if (_params.fusion_mode & MASK_USE_DRAG) {
		if (_control_status.flags.in_air
				&& !_mag_inhibit_yaw_reset_req) {
			if (!_control_status.flags.wind) {
				// reset the wind states and covariances when starting drag accel fusion
				_control_status.flags.wind = true;
				resetWindStates();
				resetWindCovariance();

			} else if (_drag_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_drag_sample_delayed)) {
				fuseDrag();

			}

		} else {
			_control_status.flags.wind = false;

		}
	}
}

void Ekf::controlMagFusion()
{
	if (_params.mag_fusion_type >= MAG_FUSE_TYPE_NONE) {

		// do not use the magnetometer and deactivate magnetic field states
		// save covariance data for re-use if currently doing 3-axis fusion
		if (_control_status.flags.mag_3D) {
			save_mag_cov_data();
			_control_status.flags.mag_3D = false;
		}
		zeroRows(P, 16, 21);
		zeroCols(P, 16, 21);
		_mag_decl_cov_reset = false;
		_control_status.flags.mag_hdg = false;

		return;
	}

	// If we are on ground, store the local position and time to use as a reference
	// Also reset the flight alignment flag so that the mag fields will be re-initialised next time we achieve flight altitude
	if (!_control_status.flags.in_air) {
		_last_on_ground_posD = _state.pos(2);
		_control_status.flags.mag_align_complete = false;
		_num_bad_flight_yaw_events = 0;
	}

	// check for new magnetometer data that has fallen behind the fusion time horizon
	// If we are using external vision data for heading then no magnetometer fusion is used
	if (!_control_status.flags.ev_yaw && !_control_status.flags.gps_yaw && _mag_data_ready) {

		// We need to reset the yaw angle after climbing away from the ground to enable
		// recovery from ground level magnetic interference.
		if (!_control_status.flags.mag_align_complete && _control_status.flags.in_air) {
			// Check if height has increased sufficiently to be away from ground magnetic anomalies
			// and request a yaw reset if not already requested.
			float terrain_vpos_estimate = isTerrainEstimateValid() ? _terrain_vpos : _last_on_ground_posD;
			_mag_yaw_reset_req |= (terrain_vpos_estimate - _state.pos(2)) > 1.5f;
		}

		// perform a yaw reset if requested by other functions
		if (_mag_yaw_reset_req && _control_status.flags.tilt_align) {
			if (!_mag_use_inhibit ) {
				if (!_control_status.flags.mag_align_complete && _control_status.flags.fixed_wing && _control_status.flags.in_air) {
					// A fixed wing vehicle can use GPS to bound yaw errors immediately after launch
					_control_status.flags.mag_align_complete = realignYawGPS();

					if (_velpos_reset_request) {
						resetToGeneralVelocity();
						resetToGeneralHorizontalPosition();
						_velpos_reset_request = false;
					}

				} else {
					_control_status.flags.mag_align_complete = resetMagHeading(_mag_sample_delayed.mag) && _control_status.flags.in_air;
				}
			}
			_control_status.flags.yaw_align = _control_status.flags.yaw_align || _control_status.flags.mag_align_complete;
			_mag_yaw_reset_req = false;
		}

		// Determine if we should use simple magnetic heading fusion which works better when there are large external disturbances
		// or the more accurate 3-axis fusion
		if (_control_status.flags.mag_fault) {
			// do no magnetometer fusion at all
			_control_status.flags.mag_hdg = false;
			_control_status.flags.mag_3D = false;

		} else if (_params.mag_fusion_type == MAG_FUSE_TYPE_AUTO || _params.mag_fusion_type == MAG_FUSE_TYPE_AUTOFW) {
			// Check if there has been enough change in horizontal velocity to make yaw observable
			// Apply hysteresis to check to avoid rapid toggling
			if (_yaw_angle_observable) {
				_yaw_angle_observable = _accel_lpf_NE.norm() > _params.mag_acc_gate;

			} else {
				_yaw_angle_observable = _accel_lpf_NE.norm() > 2.0f * _params.mag_acc_gate;
			}

			_yaw_angle_observable = _yaw_angle_observable && (_control_status.flags.gps || _control_status.flags.ev_pos); // Do we have to add ev_vel here?

			// check if there is enough yaw rotation to make the mag bias states observable
			if (!_mag_bias_observable && (fabsf(_yaw_rate_lpf_ef) > _params.mag_yaw_rate_gate)) {
				// initial yaw motion is detected
				_mag_bias_observable = true;
				_yaw_delta_ef = 0.0f;
				_time_yaw_started = _imu_sample_delayed.time_us;

			} else if (_mag_bias_observable) {
				// monitor yaw rotation in 45 deg sections.
				// a rotation of 45 deg is sufficient to make the mag bias observable
				if (fabsf(_yaw_delta_ef) > math::radians(45.0f)) {
					_time_yaw_started = _imu_sample_delayed.time_us;
					_yaw_delta_ef = 0.0f;
				}

				// require sustained yaw motion of 50% the initial yaw rate threshold
				float min_yaw_change_req =  0.5f * _params.mag_yaw_rate_gate * (1e-6f * (float)(_imu_sample_delayed.time_us - _time_yaw_started));
				_mag_bias_observable = fabsf(_yaw_delta_ef) > min_yaw_change_req;

			} else {
				_mag_bias_observable = false;
			}

			// record the last time that movement was suitable for use of 3-axis magnetometer fusion
			if (_mag_bias_observable || _yaw_angle_observable) {
				_time_last_movement = _imu_sample_delayed.time_us;
			}

			// decide whether 3-axis magnetometer fusion can be used
			bool use_3D_fusion = _control_status.flags.tilt_align && // Use of 3D fusion requires valid tilt estimates
					_control_status.flags.in_air && // don't use when on the ground because of magnetic anomalies
					_control_status.flags.mag_align_complete &&
					((_imu_sample_delayed.time_us - _time_last_movement) < 2 * 1000 * 1000); // Using 3-axis fusion for a minimum period after to allow for false negatives

			// perform switch-over
			if (use_3D_fusion) {
				if (!_control_status.flags.mag_3D) {
					// reset the mag field covariances
					zeroRows(P, 16, 21);
					zeroCols(P, 16, 21);

					// re-instate variances for the D earth axis and XYZ body axis field
					for (uint8_t index = 0; index <= 3; index ++) {
						P[index + 18][index + 18] = _saved_mag_bf_variance[index];
					}
					// re-instate the NE axis covariance sub-matrix
					for (uint8_t row = 0; row <= 1; row ++) {
						for (uint8_t col = 0; col <= 1; col ++) {
							P[row + 16][col + 16] = _saved_mag_ef_covmat[row][col];
						}
					}
				}

				// only use one type of mag fusion at the same time
				_control_status.flags.mag_3D = _control_status.flags.mag_align_complete;
				_control_status.flags.mag_hdg = !_control_status.flags.mag_3D;

			} else {
				// save covariance data for re-use if currently doing 3-axis fusion
				if (_control_status.flags.mag_3D) {
					save_mag_cov_data();
					_control_status.flags.mag_3D = false;
				}

				_control_status.flags.mag_hdg = true;
			}

			/*
			Control switch-over between only updating the mag states to updating all states
			When flying as a fixed wing aircraft, a misaligned magnetometer can cause an error in pitch/roll and accel bias estimates.
			When MAG_FUSE_TYPE_AUTOFW is selected and the vehicle is flying as a fixed wing, then magnetometer fusion is only allowed
			to access the magnetic field states.
			*/
			_control_status.flags.update_mag_states_only = (_params.mag_fusion_type == MAG_FUSE_TYPE_AUTOFW)
					&& _control_status.flags.fixed_wing;

			// For the first 5 seconds after switching to 3-axis fusion we allow the magnetic field state estimates to stabilise
			// before they are used to constrain heading drift
			_flt_mag_align_converging = ((_imu_sample_delayed.time_us - _flt_mag_align_start_time) < (uint64_t)5e6);

			if (_control_status.flags.mag_3D && _control_status_prev.flags.update_mag_states_only && !_control_status.flags.update_mag_states_only) {
				// When re-commencing use of magnetometer to correct vehicle states
				// set the field state variance to the observation variance and zero
				// the covariance terms to allow the field states re-learn rapidly
				zeroRows(P, 16, 21);
				zeroCols(P, 16, 21);
				_mag_decl_cov_reset = false;

				for (uint8_t index = 0; index <= 5; index ++) {
					P[index + 16][index + 16] = sq(_params.mag_noise);
				}

				// save covariance data for re-use when auto-switching between heading and 3-axis fusion
				save_mag_cov_data();
			}

		} else if (_params.mag_fusion_type == MAG_FUSE_TYPE_HEADING || _params.mag_fusion_type == MAG_FUSE_TYPE_INDOOR) {
			// always use heading fusion
			_control_status.flags.mag_hdg = true;

			// save covariance data for re-use if currently doing 3-axis fusion
			if (_control_status.flags.mag_3D) {
				save_mag_cov_data();
				_control_status.flags.mag_3D = false;
			}

		} else if (_params.mag_fusion_type == MAG_FUSE_TYPE_3D) {
			if (!_control_status.flags.mag_3D && _control_status.flags.yaw_align) {
				// only commence 3-axis fusion when yaw is aligned and field states set
				_control_status.flags.mag_3D = true;
			}
		} else {
			// do no magnetometer fusion at all
			_control_status.flags.mag_hdg = false;

			// save covariance data for re-use if currently doing 3-axis fusion
			if (_control_status.flags.mag_3D) {
				save_mag_cov_data();
				_control_status.flags.mag_3D = false;
			}

		}

		// if we are using 3-axis magnetometer fusion, but without external aiding, then the declination must be fused as an observation to prevent long term heading drift
		// fusing declination when gps aiding is available is optional, but recommended to prevent problem if the vehicle is static for extended periods of time
		if (_control_status.flags.mag_3D && (!_control_status.flags.gps || (_params.mag_declination_source & MASK_FUSE_DECL))) {
			_control_status.flags.mag_dec = true;

		} else {
			_control_status.flags.mag_dec = false;
		}

		// If the user has selected auto protection against indoor magnetic field errors, only use the magnetometer
		// if a yaw angle relative to true North is required for navigation. If no GPS or other earth frame aiding
		// is available, assume that we are operating indoors and the magnetometer should not be used.
		bool user_selected = (_params.mag_fusion_type == MAG_FUSE_TYPE_INDOOR);
		bool not_using_gps = !(_params.fusion_mode & MASK_USE_GPS) || !_control_status.flags.gps;
		bool not_using_evpos = !(_params.fusion_mode & MASK_USE_EVPOS) || !_control_status.flags.ev_pos;
		bool not_using_evvel = !(_params.fusion_mode & MASK_USE_EVVEL) || !_control_status.flags.ev_vel;
		bool not_selected_evyaw =  !(_params.fusion_mode & MASK_USE_EVYAW);
		if (user_selected && not_using_gps && not_using_evpos && not_using_evvel && not_selected_evyaw) {
			_mag_use_inhibit = true;
		} else {
			_mag_use_inhibit = false;
			_mag_use_not_inhibit_us = _imu_sample_delayed.time_us;
		}

		// If magnetometer use has been inhibited continuously then a yaw reset is required for a valid heading
		if (uint32_t(_imu_sample_delayed.time_us - _mag_use_not_inhibit_us) > (uint32_t)5e6) {
			_mag_inhibit_yaw_reset_req = true;
		}

		// fuse magnetometer data using the selected methods
		if (_control_status.flags.mag_3D && _control_status.flags.yaw_align) {
			if (!_mag_decl_cov_reset) {
				// After any magnetic field covariance reset event the earth field state
				// covariances need to be corrected to incorporate knowedge of the declination
				// before fusing magnetomer data to prevent rapid rotation of the earth field
				// states for the first few observations.
				fuseDeclination(0.02f);
				_mag_decl_cov_reset = true;
				fuseMag();
			} else {
				// The normal sequence is to fuse the magnetometer data first before fusing
				// declination angle at a higher uncertainty to allow some learning of
				// declination angle over time.
				fuseMag();
				if (_control_status.flags.mag_dec) {
					fuseDeclination(0.5f);
				}
			}

		} else if (_control_status.flags.mag_hdg && _control_status.flags.yaw_align) {
			// fusion of an Euler yaw angle from either a 321 or 312 rotation sequence
			fuseHeading();

		} else {
			// do no fusion at all
		}
	}
}

void Ekf::controlFakePosFusion()
{
	// if we aren't doing any aiding, fake position measurements at the last known position to constrain drift
	// Coincide fake measurements with baro data for efficiency with a minimum fusion rate of 5Hz

	if (!_control_status.flags.gps &&
	    !_control_status.flags.opt_flow &&
	    !_control_status.flags.ev_pos &&
	    !_control_status.flags.ev_vel &&
	    !(_control_status.flags.fuse_aspd && _control_status.flags.fuse_beta)) {

		// We now need to use a synthetic position observation to prevent unconstrained drift of the INS states.
		_using_synthetic_position = true;

		// Fuse synthetic position observations every 200msec
		if (((_time_last_imu - _time_last_fake_pos) > (uint64_t)2e5) || _fuse_height) {

			Vector3f fake_pos_obs_var{};
			Vector2f fake_pos_innov_gate{};


			// Reset position and velocity states if we re-commence this aiding method
			if ((_time_last_imu - _time_last_fake_pos) > (uint64_t)4e5) {
				resetToGeneralHorizontalPosition();
				resetToGeneralVelocity();
				_fuse_hpos_as_odom = false;

				if (_time_last_fake_pos != 0) {
					ECL_WARN_TIMESTAMPED("stopping navigation");
				}

			}
			_time_last_fake_pos = _time_last_imu;

			if (_control_status.flags.in_air && _control_status.flags.tilt_align) {
				fake_pos_obs_var(0) = fake_pos_obs_var(1) = fmaxf(_params.pos_noaid_noise, _params.gps_pos_noise);

			} else {
				fake_pos_obs_var(0) = fake_pos_obs_var(1) = 0.5f;
			}

			_gps_pos_innov(0) = _state.pos(0) - _last_known_posNE(0);
			_gps_pos_innov(1) = _state.pos(1) - _last_known_posNE(1);

			// glitch protection is not required so set gate to a large value
			fake_pos_innov_gate(0) = 100.0f;

			fuseHorizontalPosition(_gps_pos_innov, fake_pos_innov_gate, fake_pos_obs_var,
						_gps_pos_innov_var, _gps_pos_test_ratio);
		}

	} else {
		_using_synthetic_position = false;
	}

}

void Ekf::controlAuxVelFusion()
{
	if(_aux_vel_data_ready)
	{
		if ((_params.fusion_mode & MASK_USE_AUXVEL) && !_control_status.flags.aux_vel
			&& (_time_last_imu - _time_last_auxvel) < (2 * EV_MAX_INTERVAL)){
			resetToAuxiliarHorizontalVelocity();
			startAuxVelFusion();
		}

		if (_control_status.flags.aux_vel) {

			Vector2f aux_vel_innov_gate{};
			Vector3f aux_vel_obs_var{};

			// TODO: Add vertical velocity
			_aux_vel_innov(0) = _state.vel(0) - _auxvel_sample_delayed.vel(0);
			_aux_vel_innov(1) = _state.vel(1) - _auxvel_sample_delayed.vel(1);
			aux_vel_innov_gate(0) = _params.auxvel_gate;
			aux_vel_obs_var(0) = _auxvel_sample_delayed.velVar(0);
			aux_vel_obs_var(1) = _auxvel_sample_delayed.velVar(1);

			fuseHorizontalVelocity(_aux_vel_innov, aux_vel_innov_gate, aux_vel_obs_var,
					_aux_vel_innov_var, _aux_vel_test_ratio);

		}
	}
	else if(_control_status.flags.aux_vel && isTimedOut(_time_last_auxvel, 1e6f)){
		stopAuxVelFusion();
		ECL_INFO_TIMESTAMPED("stop auxiliar velocity fusion - data timeout");
	}
}

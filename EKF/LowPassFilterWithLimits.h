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
#pragma once

#include <mathlib/mathlib.h>
#include <ecl.h>

/**
 * Lowpass filter with limits on filter state and filter input
 * @author Kamil Ritz <ka.ritz@hotmail.com>
 */

namespace estimator {

class LowPassFilterWithLimits
{
public:
	// LowPassFilterWithLimits(float max_input_dev = NAN; float min_value = NAN; float max_value = NAN):
	// 	_max_input_deviation(max_input_dev), _min_value(min_value), _max_value(max_value)
	// 	{reset();};
	LowPassFilterWithLimits() = default;
	~LowPassFilterWithLimits() = default;

	void reset(float val = 0.f, uint64_t t_us = 0)
	{
		_x = val;
		_time_last_update_us = t_us;
	}

	/**
	 * Update the filter with a new value and returns the filtered state
	 * The input value is constrained to be close to current state by the limit sets in setSpikeLimit
	 * @param val new input
	 * @param alpha normalized weight of the new input
	 * @return filtered output
	 */
	float update(float val, float alpha)
	{
		float input_constrained;

		if(ISFINITE(_max_input_deviation)){
			input_constrained = math::constrain(val, _x - _max_input_deviation, _x + _max_input_deviation);
		}else{
			input_constrained = val;
		}
		float beta = 1.f - alpha;

		_x = beta * _x + alpha * input_constrained;

		if (ISFINITE(_min_value))
		{
			_x = math::max(_min_value, _x);
		}
		if (ISFINITE(_max_value))
		{
			_x = math::min(_max_value, _x);
		}

		return _x;
	}

	float getValue()
	{
		return _x;
	}

	float getMaxValue()
	{
		return _max_value;
	}

	float getMinValue()
	{
		return _min_value;
	}

	float getMaxInputDeviation()
	{
		return _max_input_deviation;
	}

	uint64_t getLastUpdateTime()
	{
		return _time_last_update_us;
	}

	/**
	 * Helper function to compute alpha from dt and time constant tau
	 * Under assumption that sampling time is much smaller than time constant
	 * @param dt sampling time in microseconds
	 * @param tau_inv inverse of the time constant of the filter (1/s)
	 * @return alpha, the normalized weight of a new measurement
	 */
	static float computeAlphaFromDtAndTauInv(float dt, float tau_inv)
	{
		return math::constrain(dt * tau_inv * 1e-6f, 0.0f, 1.0f);
	}

	/**
	 * Helper function to compute alpha from current timestamp and time constant tau
	 * Under assumption that sampling time is much smaller than time constant
	 * @param t_curr_us timestamp of current update (us)
	 * @param tau_inv inverse of the time constant of the filter (1/s)
	 * @return alpha, the normalized weight of a new measurement
	 */
	float computeAlphaFromTimestampAndTauInv(uint64_t t_curr_us, float tau_inv)
	{
		float alpha = math::constrain(float(t_curr_us - _time_last_update_us) * 1e-6f * tau_inv,0.0f,1.0f);
		_time_last_update_us = t_curr_us;
		return alpha;
	}

	void setMaxInputDeviation(float max_input_deviation)
	{
		if(max_input_deviation>=0.0f)
			_max_input_deviation = max_input_deviation;
		else{
			_max_input_deviation = 1.0f;
		}
	}

	void setStateLimits(float min_value, float max_value)
	{
		if(min_value <= max_value){
			_min_value = min_value;
			_max_value = max_value;
		}else{
			_min_value = 1.0f;
			_max_value = 1.0f;
		}

	}



private:
	float _x{0.0f}; ///< current state of the filter
	uint64_t _time_last_update_us{0}; ///< timestamp of last filter update (us)
	float _max_input_deviation{NAN}; ///< max deviation of the input to the state of the filter, if = NAN no limits apply
	float _min_value{NAN}; ///< min value the state can take, if = NAN no limits apply
	float _max_value{NAN}; ///< max value the state can take, if = NAN no limits apply
};

}

/****************************************************************************
 *
 *   Copyright (c) 2019 ECL Development Team. All rights reserved.
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

#include <gtest/gtest.h>
#include <cmath>
#include "EKF/LowPassFilterWithLimits.h"


class EkfLowPassFilterWithLimitsTest : public ::testing::Test {
 public:

	estimator::LowPassFilterWithLimits _filter{};

	void SetUp() override
	{

	}

	void TearDown() override
	{

	}
};

TEST_F(EkfLowPassFilterWithLimitsTest, testConstructor)
{
	// GIVEN: a created filter object

	// THEN: state should be initialized to zero
	EXPECT_FLOAT_EQ(0.0f, _filter.getValue());

	// THEN: time of last update should be initialised to zero
	EXPECT_EQ(0, _filter.getLastUpdateTime());

	// THEN: Limits should be set to NAN
	EXPECT_EQ(false, ISFINITE(_filter.getMaxInputDeviation()));
	EXPECT_EQ(false, ISFINITE(_filter.getMaxValue()));
	EXPECT_EQ(false, ISFINITE(_filter.getMinValue()));

	// WHEN: computing alpha after 1s
	float tau_inv = 1.0f/555.55f;
	float alpha = _filter.computeAlphaFromTimestampAndTauInv(1e6, tau_inv);
	// THEN: alpha should be equal to the inverse of the timeconstant
	EXPECT_FLOAT_EQ(tau_inv, alpha);
}

TEST_F(EkfLowPassFilterWithLimitsTest, testFilterConstantCalculation)
{
	// GIVEN: a reset filter object
	_filter.reset(0.0f,(uint64_t)1e6);

	// WHEN: computing alpha with the last timestamp
	float alpha = _filter.computeAlphaFromTimestampAndTauInv((uint64_t)1e6, 0.5f);

	// THEN: alpha should be zero
	EXPECT_FLOAT_EQ(0.0f, alpha);
}

TEST_F(EkfLowPassFilterWithLimitsTest, testFilterUpdate)
{
	// GIVEN: a fix alpha value
	float alpha = 0.3f;
	float beta = 1.0f - alpha;
	// WHEN: applying a step to filter without a spike limit
	_filter.reset(0);
	for (int i = 0; i < 5; i++)
	{
		_filter.update(1.0f, alpha);
	}

	// THEN: alpha should be equal to precomputed value
	float exp_filter_state = beta * beta * beta *beta + beta * beta * beta + beta * beta + beta + 1;
	exp_filter_state *= alpha;
	EXPECT_FLOAT_EQ(exp_filter_state, _filter.getValue());
}

TEST_F(EkfLowPassFilterWithLimitsTest, testFilterConvergence)
{
	// GIVEN: a fix alpha value
	float alpha = 0.3f;
	// WHEN: applying a step to filter for long time even with harsh spike limit
	_filter.reset(0);
	for (int i = 0; i < 1000; i++)
	{
		_filter.update(1.0f, alpha);
	}

	// THEN: the filter value should converge to the input
	EXPECT_NEAR(1.0f, _filter.getValue(), 1e-5);

	// WHEN: applying a step to filter for long time even with harsh spike limit
	_filter.reset(0);
	for (int i = 0; i < 1000; i++)
	{
		_filter.update(-1.0f, alpha);
	}

	// THEN: the filter value should converge to the input
	EXPECT_NEAR(-1.0f, _filter.getValue(), 1e-5);
}

TEST_F(EkfLowPassFilterWithLimitsTest, testStateLimits)
{
	// GIVEN: a zero filter state

	// WHEN: applying a big positive input that is greater than the max_value state limit
	float alpha = 0.5f;
	float max_value = 1.0f;
	float min_value = -1.0f;
	_filter.reset(0);
	_filter.setStateLimits(min_value,max_value);

	for (int i = 0; i < 1000; i++)
	{
		_filter.update(2.0f*max_value, alpha);
	}

	// THEN: the filter value should be the max_value state limit
	EXPECT_FLOAT_EQ(max_value, _filter.getValue());

	// WHEN: applying a big negative input that is smaller than the min_value state limit
	for (int i = 0; i < 1000; i++)
	{
		_filter.update(2.0f*min_value, alpha);
	}

	// THEN: the filter value should be the max_value state limit
	EXPECT_FLOAT_EQ(min_value, _filter.getValue());
}

TEST_F(EkfLowPassFilterWithLimitsTest, testInputLimits)
{
	// GIVEN: a fix alpha value and zero filter state
	float alpha = 1.0f;
	float limit = 0.5f;
	_filter.reset(0);

	// WHEN: applying an input that deviates to much from the current state
	_filter.setMaxInputDeviation(limit);
	_filter.update(1.0f, alpha);

	// THEN: the filter input should be limited
	EXPECT_FLOAT_EQ(limit, _filter.getValue());
}

#pragma once
// numeric_limits.hpp: definition of numeric_limits for the balanced tau number system
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

namespace std {

template <unsigned nbits, unsigned rbits, typename bt>
class numeric_limits< sw::universal::tau<nbits, rbits, bt> > {
public:
	using TAU = sw::universal::tau<nbits, rbits, bt>;
	static constexpr bool is_specialized = true;
	static constexpr TAU min() { // return minimum positive value
		return TAU(sw::universal::SpecificValue::minpos);
	}
	static constexpr TAU max() { // return maximum value
		return TAU(sw::universal::SpecificValue::maxpos);
	}
	static constexpr TAU lowest() { // return most negative value
		return TAU(sw::universal::SpecificValue::maxneg);
	}
	static constexpr TAU epsilon() { // smallest increment from 1.0 (uniform resolution grid)
		return TAU(sw::universal::SpecificValue::minpos);
	}
	static constexpr TAU round_error() { // return largest rounding error
		return TAU(0.5);
	}
	static constexpr TAU denorm_min() {  // return minimum denormalized value
		return TAU(sw::universal::SpecificValue::minpos);
	}
	static constexpr TAU infinity() { // no infinity encoding: saturate
		return TAU(sw::universal::SpecificValue::maxpos);
	}
	static constexpr TAU quiet_NaN()     { return TAU(sw::universal::SpecificValue::zero); }
	static constexpr TAU signaling_NaN() { return TAU(sw::universal::SpecificValue::zero); }

	static constexpr int digits       = static_cast<int>(nbits);
	static constexpr int digits10     = static_cast<int>(nbits / 2.0902f); // log10(tau) ~= 0.4178 -> 1/0.4178
	static constexpr int max_digits10 = digits10 + 1;
	static constexpr bool is_signed   = true;
	static constexpr bool is_integer  = false;
	static constexpr bool is_exact    = false;
	static constexpr int radix        = 2; // irrational base tau ~= 2.618 has no integer radix

	static constexpr int min_exponent   = -static_cast<int>(rbits);
	static constexpr int min_exponent10 = static_cast<int>(min_exponent * 0.4178f);
	static constexpr int max_exponent   = static_cast<int>(nbits) - 1 - static_cast<int>(rbits);
	static constexpr int max_exponent10 = static_cast<int>(max_exponent * 0.4178f);
	static constexpr bool has_infinity  = false;
	static constexpr bool has_quiet_NaN = false;
	static constexpr bool has_signaling_NaN = false;
	static constexpr float_denorm_style has_denorm = denorm_absent;
	static constexpr bool has_denorm_loss = false;

	static constexpr bool is_iec559  = false;
	static constexpr bool is_bounded = true;
	static constexpr bool is_modulo  = false;
	static constexpr bool traps      = false;
	static constexpr bool tinyness_before = false;
	static constexpr float_round_style round_style = round_to_nearest;
};

}

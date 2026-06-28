#pragma once
// mathlib.hpp: definition of mathematical functions for the balanced tau number system
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <cmath>

// The balanced tau system uses an irrational base, so the elementary functions
// are evaluated in native floating-point and the result is rounded back into the
// nearest representable tau value. This gives a complete "native type" math
// surface without claiming bit-exact transcendental results.

namespace sw { namespace universal {

	template<unsigned nbits, unsigned rbits, typename bt>
	using TauT = tau<nbits, rbits, bt>;

	// classification
	template<unsigned nbits, unsigned rbits, typename bt> inline bool isnan(const tau<nbits, rbits, bt>&)    { return false; }
	template<unsigned nbits, unsigned rbits, typename bt> inline bool isinf(const tau<nbits, rbits, bt>&)    { return false; }
	template<unsigned nbits, unsigned rbits, typename bt> inline bool isfinite(const tau<nbits, rbits, bt>&) { return true; }
	template<unsigned nbits, unsigned rbits, typename bt> inline bool isnormal(const tau<nbits, rbits, bt>& v) { return !v.iszero(); }

	// abs is declared in tau_impl.hpp; provide fabs as an alias
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> fabs(const tau<nbits, rbits, bt>& v) { return abs(v); }

	// copysign returns a value with the magnitude of a and the sign of b
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> copysign(const tau<nbits, rbits, bt>& a, const tau<nbits, rbits, bt>& b) {
		tau<nbits, rbits, bt> mag = abs(a);
		return b.isneg() ? -mag : mag;
	}

	// square root
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> sqrt(const tau<nbits, rbits, bt>& v) {
		if (v.isneg()) {
#if TAU_THROW_ARITHMETIC_EXCEPTION
			throw tau_negative_sqrt_arg();
#else
			return tau<nbits, rbits, bt>(0);
#endif
		}
		return tau<nbits, rbits, bt>(std::sqrt(v.to_double()));
	}

#define TAU_UNARY_MATH_FUNCTION(NAME)                                                       \
	template<unsigned nbits, unsigned rbits, typename bt>                                    \
	inline tau<nbits, rbits, bt> NAME(const tau<nbits, rbits, bt>& v) {                      \
		return tau<nbits, rbits, bt>(std::NAME(v.to_double()));                              \
	}

	TAU_UNARY_MATH_FUNCTION(exp)
	TAU_UNARY_MATH_FUNCTION(exp2)
	TAU_UNARY_MATH_FUNCTION(log)
	TAU_UNARY_MATH_FUNCTION(log2)
	TAU_UNARY_MATH_FUNCTION(log10)
	TAU_UNARY_MATH_FUNCTION(sin)
	TAU_UNARY_MATH_FUNCTION(cos)
	TAU_UNARY_MATH_FUNCTION(tan)
	TAU_UNARY_MATH_FUNCTION(asin)
	TAU_UNARY_MATH_FUNCTION(acos)
	TAU_UNARY_MATH_FUNCTION(atan)
	TAU_UNARY_MATH_FUNCTION(sinh)
	TAU_UNARY_MATH_FUNCTION(cosh)
	TAU_UNARY_MATH_FUNCTION(tanh)
	TAU_UNARY_MATH_FUNCTION(floor)
	TAU_UNARY_MATH_FUNCTION(ceil)
	TAU_UNARY_MATH_FUNCTION(trunc)
	TAU_UNARY_MATH_FUNCTION(round)

#undef TAU_UNARY_MATH_FUNCTION

	// power and two-argument functions
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> pow(const tau<nbits, rbits, bt>& a, const tau<nbits, rbits, bt>& b) {
		return tau<nbits, rbits, bt>(std::pow(a.to_double(), b.to_double()));
	}
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> atan2(const tau<nbits, rbits, bt>& y, const tau<nbits, rbits, bt>& x) {
		return tau<nbits, rbits, bt>(std::atan2(y.to_double(), x.to_double()));
	}
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> hypot(const tau<nbits, rbits, bt>& x, const tau<nbits, rbits, bt>& y) {
		return tau<nbits, rbits, bt>(std::hypot(x.to_double(), y.to_double()));
	}
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> fmod(const tau<nbits, rbits, bt>& x, const tau<nbits, rbits, bt>& y) {
		return tau<nbits, rbits, bt>(std::fmod(x.to_double(), y.to_double()));
	}

	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> min(const tau<nbits, rbits, bt>& a, const tau<nbits, rbits, bt>& b) { return (a < b) ? a : b; }
	template<unsigned nbits, unsigned rbits, typename bt>
	inline tau<nbits, rbits, bt> max(const tau<nbits, rbits, bt>& a, const tau<nbits, rbits, bt>& b) { return (a > b) ? a : b; }

}} // namespace sw::universal

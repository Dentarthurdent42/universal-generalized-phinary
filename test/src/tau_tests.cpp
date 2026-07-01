// tau_tests.cpp: regression tests for the balanced tau (phi^2) number system
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <cstdint>
#include <limits>
#include <string>

#include <universal/number/tau/tau.hpp>

using namespace sw::universal;

// The base tau = phi^2 = (3 + sqrt(5))/2
static constexpr double TAU = 2.6180339887498948482;

// Integers round-trip exactly within the trit budget.
TEST(TauTest, IntegerRoundTrip) {
	using Tau = tau<20, 8>;
	for (int n = -50; n <= 50; ++n) {
		Tau v = n;
		EXPECT_NEAR(v.to_double(), double(n), 1e-9) << "n = " << n;
	}
}

// The defining identity of the base: tau + 1/tau = 3 (from tau^2 - 3 tau + 1 = 0).
TEST(TauTest, DefiningIdentity) {
	using Tau = tau<24, 12>;
	Tau t = TAU;
	Tau recip = 1.0 / TAU;
	Tau three = t + recip;
	EXPECT_NEAR(three.to_double(), 3.0, 1e-3);
}

// Balanced digit set yields a symmetric range about zero.
TEST(TauTest, SymmetricRange) {
	using Tau = tau<16, 6>;
	Tau maxpos(SpecificValue::maxpos);
	Tau maxneg(SpecificValue::maxneg);
	EXPECT_EQ(maxneg, -maxpos);
	EXPECT_DOUBLE_EQ(maxneg.to_double(), -maxpos.to_double());
}

// Negation is exact (trit-by-trit), so -(-x) == x.
TEST(TauTest, Negation) {
	using Tau = tau<20, 8>;
	for (int n = -20; n <= 20; ++n) {
		Tau v = n;
		EXPECT_EQ(-(-v), v) << "n = " << n;
		EXPECT_NEAR((-v).to_double(), -double(n), 1e-9) << "n = " << n;
	}
}

// Basic arithmetic on integer operands stays exact.
TEST(TauTest, IntegerArithmetic) {
	using Tau = tau<24, 8>;
	Tau a = 6, b = 7;
	EXPECT_NEAR((a + b).to_double(), 13.0, 1e-9);
	EXPECT_NEAR((a - b).to_double(), -1.0, 1e-9);
	EXPECT_NEAR((a * b).to_double(), 42.0, 1e-9);
	EXPECT_NEAR((Tau(42) / Tau(6)).to_double(), 7.0, 1e-9);
}

// Quadratic formula: x^2 - 5x + 6 = 0 has roots 2 and 3.
TEST(TauTest, SolveSimpleQuadratic) {
	using Tau = tau<32, 12>;
	Tau a(1), b(-5), c(6);
	Tau discriminant = b * b - Tau(4) * a * c;
	Tau sqrt_disc = sqrt(discriminant);
	Tau root1 = (-b - sqrt_disc) / (Tau(2) * a);
	Tau root2 = (-b + sqrt_disc) / (Tau(2) * a);
	EXPECT_NEAR(root1.to_double(), 2.0, 1e-3);
	EXPECT_NEAR(root2.to_double(), 3.0, 1e-3);
}

// Ordering reflects the represented real value.
TEST(TauTest, Comparison) {
	using Tau = tau<20, 8>;
	EXPECT_TRUE(Tau(-3) < Tau(2));
	EXPECT_TRUE(Tau(5) > Tau(4));
	EXPECT_TRUE(Tau(7) >= Tau(7));
	EXPECT_TRUE(Tau(0) <= Tau(1));
	EXPECT_TRUE(Tau(2) == Tau(2));
	EXPECT_TRUE(Tau(2) != Tau(3));
}

// Selectors: sign, zero, scale.
TEST(TauTest, Selectors) {
	using Tau = tau<16, 6>;
	EXPECT_TRUE(Tau(0).iszero());
	EXPECT_FALSE(Tau(1).iszero());
	EXPECT_TRUE(Tau(-4).isneg());
	EXPECT_FALSE(Tau(4).isneg());
	EXPECT_TRUE(Tau(4).ispos());
	// scale is the exponent of the most significant non-zero trit
	EXPECT_EQ(Tau(1).scale(), 0);
}

// Values can be used in constexpr expressions (C++23 requirement for arithmetic types).
TEST(TauTest, ConstexprUsage) {
	using Tau = tau<16, 6>;
	constexpr Tau one(1);
	static_assert(one.nbits == 16, "nbits mismatch");
	static_assert(one.rbits == 6, "rbits mismatch");
	constexpr long double v = one.to_long_double();
	static_assert(v > 0.99L && v < 1.01L, "constexpr value of one");

	constexpr Tau two(2);
	constexpr Tau sum = one + two;
	static_assert(sum.to_long_double() > 2.99L && sum.to_long_double() < 3.01L, "constexpr 1+2");
	SUCCEED();
}

// numeric_limits specialization is wired up.
TEST(TauTest, NumericLimits) {
	using Tau = tau<16, 6>;
	EXPECT_TRUE(std::numeric_limits<Tau>::is_specialized);
	EXPECT_TRUE(std::numeric_limits<Tau>::is_signed);
	EXPECT_GT(std::numeric_limits<Tau>::max().to_double(), 0.0);
	EXPECT_LT(std::numeric_limits<Tau>::lowest().to_double(), 0.0);
	EXPECT_GT(std::numeric_limits<Tau>::epsilon().to_double(), 0.0);
}

// type traits identify tau types.
TEST(TauTest, Traits) {
	EXPECT_TRUE((is_tau< tau<16, 6, uint32_t> >));
	EXPECT_FALSE(is_tau<double>);
}

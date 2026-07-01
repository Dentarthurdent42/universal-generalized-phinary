#pragma once
// tau_impl.hpp: implementation of the balanced tau (phi^2) number system
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>

// supporting types and functions
#include <universal/native/ieee754.hpp>
#include <universal/number/shared/specific_value_encoding.hpp>
#include <universal/internal/blockbinary/blockbinary.hpp>

// the exception types
#include <universal/number/tau/exceptions.hpp>

/*
The balanced tau number system is the "ternary phi-representation" of
Agaian and Zhou's Generalized Phi Number System (GPNS): a positional number
system whose base is

    tau = phi^2 = (3 + sqrt(5)) / 2 ~= 2.6180339887...

and whose digits are the balanced-ternary set { -1, 0, +1 }.  A value is

    X = sum_p  d_p * tau^(p - rbits),   d_p in {-1,0,+1}

where p in [0, nbits) indexes the trits, the rbits least-significant trits
sit to the right of the radix point, and the remaining (nbits - rbits) trits
form the integer part.

Because the digit set is symmetric, the system represents both positive and
negative numbers without a dedicated sign field, and the representable range
is symmetric about zero (maxneg() == -maxpos()).

tau is an irrational base, so the set of representable values is not closed
under the field operations.  Arithmetic is therefore performed in extended
precision and the result is rounded to the nearest representable tau value
(saturating at the extremes of the range).  This mirrors how every finite
number system rounds; the trit encoding remains the canonical source of truth.
*/

namespace sw { namespace universal {

// enable/disable throwing arithmetic exceptions
#if !defined(TAU_THROW_ARITHMETIC_EXCEPTION)
#define TAU_THROW_ARITHMETIC_EXCEPTION 0
#endif

// trit-packing helpers store each balanced trit in two bits: bit 2p is the
// magnitude (1 == non-zero) and bit 2p+1 is the sign (1 == negative).

// template class representing a value in the balanced base-tau number system
template<unsigned _nbits, unsigned _rbits, typename bt = std::uint32_t>
class tau {
public:
	static_assert(_nbits >= _rbits, "tau configuration error: nbits must be greater or equal to rbits");
	static_assert(_nbits >= 1, "tau configuration error: nbits must be at least 1");

	static constexpr unsigned nbits = _nbits;   // total number of trits
	static constexpr unsigned rbits = _rbits;   // number of fractional trits
	static constexpr unsigned fbits = _rbits;   // symmetry with other Universal types
	static constexpr unsigned ibits = _nbits - _rbits; // number of integer trits
	typedef bt BlockType;

	// the irrational base: tau = phi^2 = (3 + sqrt(5))/2, and the golden ratio phi
	static constexpr long double tau_base = 2.6180339887498948482045868343656381177L;
	static constexpr long double phi      = 1.6180339887498948482045868343656381177L;

	// integer power of the base tau, usable in constexpr context
	static constexpr long double ipow(int e) noexcept {
		long double r = 1.0L;
		if (e >= 0) { for (int i = 0; i < e; ++i) r *= tau_base; }
		else        { for (int i = 0; i < -e; ++i) r /= tau_base; }
		return r;
	}
	// largest representable magnitude: all trits set to +1
	static constexpr long double maxpos_value() noexcept {
		long double v = 0.0L;
		for (unsigned p = 0; p < nbits; ++p) v += ipow(static_cast<int>(p) - static_cast<int>(rbits));
		return v;
	}
	// smallest representable positive magnitude: lowest trit set to +1
	static constexpr long double minpos_value() noexcept {
		return ipow(-static_cast<int>(rbits));
	}

	// constructors
	constexpr tau() noexcept { _block.clear(); }

	constexpr tau(const tau&) noexcept = default;
	constexpr tau(tau&&) noexcept = default;
	constexpr tau& operator=(const tau&) noexcept = default;
	constexpr tau& operator=(tau&&) noexcept = default;

	// converting constructors from native types
	constexpr tau(signed char  iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(short        iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(int          iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(long         iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(long long    iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(unsigned int       iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(unsigned long      iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(unsigned long long iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(float       iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(double      iv) noexcept { _block.clear(); convert(static_cast<long double>(iv)); }
	constexpr tau(long double iv) noexcept { _block.clear(); convert(iv); }

	// specific value constructor
	constexpr tau(const SpecificValue code) noexcept {
		_block.clear();
		switch (code) {
		case SpecificValue::infpos:
		case SpecificValue::maxpos:
			maxpos();
			break;
		case SpecificValue::minpos:
			minpos();
			break;
		case SpecificValue::minneg:
			minneg();
			break;
		case SpecificValue::infneg:
		case SpecificValue::maxneg:
			maxneg();
			break;
		case SpecificValue::qnan:
		case SpecificValue::snan:
		case SpecificValue::nar:
		case SpecificValue::zero:
		default:
			setzero();
			break;
		}
	}

	// assignment operators from native types
	constexpr tau& operator=(signed char  rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(short        rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(int          rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(long         rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(long long    rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(unsigned int       rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(unsigned long      rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(unsigned long long rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(float       rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(double      rhs) noexcept { return convert(static_cast<long double>(rhs)); }
	constexpr tau& operator=(long double rhs) noexcept { return convert(rhs); }

	// arithmetic operators
	// prefix negation: balanced digits negate trit-by-trit, which is exact
	constexpr tau operator-() const noexcept {
		tau negated;
		for (unsigned p = 0; p < nbits; ++p) negated.setTrit(p, -trit(p));
		return negated;
	}

	// in-place arithmetic assignment operators (extended precision, then round)
	constexpr tau& operator+=(const tau& rhs) noexcept { return convert(to_long_double() + rhs.to_long_double()); }
	constexpr tau& operator+=(double rhs)      noexcept { return *this += tau(rhs); }
	constexpr tau& operator-=(const tau& rhs) noexcept { return convert(to_long_double() - rhs.to_long_double()); }
	constexpr tau& operator-=(double rhs)      noexcept { return *this -= tau(rhs); }
	constexpr tau& operator*=(const tau& rhs) noexcept { return convert(to_long_double() * rhs.to_long_double()); }
	constexpr tau& operator*=(double rhs)      noexcept { return *this *= tau(rhs); }
	constexpr tau& operator/=(const tau& rhs) noexcept {
		if (rhs.iszero()) {
#if TAU_THROW_ARITHMETIC_EXCEPTION
			throw tau_divide_by_zero();
#else
			// saturate to the extreme of the range that matches the dividend's sign
			bool neg = isneg();
			setzero();
			if (neg) maxneg(); else maxpos();
			return *this;
#endif
		}
		return convert(to_long_double() / rhs.to_long_double());
	}
	constexpr tau& operator/=(double rhs) noexcept { return *this /= tau(rhs); }

	// prefix/postfix operators (unit increment/decrement)
	constexpr tau& operator++()    noexcept { return convert(to_long_double() + 1.0L); }
	constexpr tau  operator++(int) noexcept { tau tmp(*this); operator++(); return tmp; }
	constexpr tau& operator--()    noexcept { return convert(to_long_double() - 1.0L); }
	constexpr tau  operator--(int) noexcept { tau tmp(*this); operator--(); return tmp; }

	// modifiers
	constexpr void clear()   noexcept { _block.clear(); }
	constexpr void setzero() noexcept { _block.clear(); }
	constexpr tau& maxpos() noexcept { for (unsigned p = 0; p < nbits; ++p) setTrit(p,  1); return *this; }
	constexpr tau& maxneg() noexcept { for (unsigned p = 0; p < nbits; ++p) setTrit(p, -1); return *this; }
	constexpr tau& minpos() noexcept { _block.clear(); setTrit(0,  1); return *this; }
	constexpr tau& minneg() noexcept { _block.clear(); setTrit(0, -1); return *this; }

	// set/get a single balanced trit at position p in [0, nbits)
	constexpr void setTrit(unsigned p, int d) noexcept {
		_block.setbit(2 * p,     d != 0);
		_block.setbit(2 * p + 1, d <  0);
	}
	constexpr int trit(unsigned p) const noexcept {
		if (!_block.at(2 * p)) return 0;
		return _block.at(2 * p + 1) ? -1 : 1;
	}

	// selectors
	constexpr bool iszero() const noexcept {
		for (unsigned p = 0; p < nbits; ++p) if (trit(p) != 0) return false;
		return true;
	}
	constexpr bool isneg() const noexcept { return sign(); }
	constexpr bool ispos() const noexcept { return !iszero() && !sign(); }
	constexpr bool isinf() const noexcept { return false; }
	constexpr bool isnan() const noexcept { return false; }
	// the sign is the sign of the most-significant non-zero trit
	constexpr bool sign() const noexcept {
		for (int p = static_cast<int>(nbits) - 1; p >= 0; --p) {
			int d = trit(static_cast<unsigned>(p));
			if (d != 0) return d < 0;
		}
		return false;
	}
	// exact base-tau scale: exponent of the most-significant non-zero trit
	constexpr int scale() const noexcept {
		for (int p = static_cast<int>(nbits) - 1; p >= 0; --p) {
			if (trit(static_cast<unsigned>(p)) != 0) return p - static_cast<int>(rbits);
		}
		return 0;
	}

	constexpr long double maxpos_bound() const noexcept { return maxpos_value(); }

	// access to the underlying storage
	constexpr const blockbinary<2 * nbits, bt>& bits() const noexcept { return _block; }

	// value extraction
	constexpr long double to_long_double() const noexcept {
		long double v = 0.0L;
		for (unsigned p = 0; p < nbits; ++p) {
			int d = trit(p);
			if (d != 0) v += static_cast<long double>(d) * ipow(static_cast<int>(p) - static_cast<int>(rbits));
		}
		return v;
	}
	constexpr double to_double() const noexcept { return static_cast<double>(to_long_double()); }
	constexpr float  to_float()  const noexcept { return static_cast<float>(to_long_double()); }

	explicit constexpr operator long double() const noexcept { return to_long_double(); }
	explicit constexpr operator double()      const noexcept { return to_double(); }
	explicit constexpr operator float()       const noexcept { return to_float(); }
	explicit constexpr operator int()         const noexcept { return static_cast<int>(to_long_double()); }
	explicit constexpr operator long long()   const noexcept { return static_cast<long long>(to_long_double()); }

private:
	blockbinary<2 * nbits, bt> _block;

	// encode a real value into the nearest representable balanced base-tau number.
	// Greedy nearest-digit conversion from the most- to the least-significant trit,
	// saturating to the representable range.
	constexpr tau& convert(long double x) noexcept {
		_block.clear();
		if (x == 0.0L) return *this;
		const long double maxv = maxpos_value();
		if (x >  maxv) x =  maxv;
		if (x < -maxv) x = -maxv;
		long double r = x;
		for (int p = static_cast<int>(nbits) - 1; p >= 0; --p) {
			long double w = ipow(p - static_cast<int>(rbits));
			long double q = r / w;
			long long n = (q >= 0.0L) ? static_cast<long long>(q + 0.5L)
			                          : static_cast<long long>(q - 0.5L);
			int d = (n > 0) ? 1 : (n < 0 ? -1 : 0);
			if (d != 0) {
				setTrit(static_cast<unsigned>(p), d);
				r -= static_cast<long double>(d) * w;
			}
		}
		return *this;
	}

	// befriend the operators that need access to the private encoding
	template<unsigned nnbits, unsigned nrbits, typename nbt>
	friend std::ostream& operator<<(std::ostream&, const tau<nnbits, nrbits, nbt>&);
	template<unsigned nnbits, unsigned nrbits, typename nbt>
	friend std::istream& operator>>(std::istream&, tau<nnbits, nrbits, nbt>&);
};

////////////////////// stream operators
template<unsigned nnbits, unsigned nrbits, typename nbt>
inline std::ostream& operator<<(std::ostream& ostr, const tau<nnbits, nrbits, nbt>& v) {
	return ostr << v.to_double();
}
template<unsigned nnbits, unsigned nrbits, typename nbt>
inline std::istream& operator>>(std::istream& istr, tau<nnbits, nrbits, nbt>& v) {
	double d{ 0 };
	istr >> d;
	v = d;
	return istr;
}

////////////////////// comparison operators (value based)
template<unsigned nbits, unsigned rbits, typename bt>
constexpr bool operator==(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) {
	return lhs.to_long_double() == rhs.to_long_double();
}
template<unsigned nbits, unsigned rbits, typename bt>
constexpr bool operator!=(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) { return !operator==(lhs, rhs); }
template<unsigned nbits, unsigned rbits, typename bt>
constexpr bool operator< (const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) {
	return lhs.to_long_double() < rhs.to_long_double();
}
template<unsigned nbits, unsigned rbits, typename bt>
constexpr bool operator> (const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) { return  operator< (rhs, lhs); }
template<unsigned nbits, unsigned rbits, typename bt>
constexpr bool operator<=(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) { return !operator> (lhs, rhs); }
template<unsigned nbits, unsigned rbits, typename bt>
constexpr bool operator>=(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) { return !operator< (lhs, rhs); }

////////////////////// binary arithmetic operators
template<unsigned nbits, unsigned rbits, typename bt>
constexpr tau<nbits, rbits, bt> operator+(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) {
	tau<nbits, rbits, bt> sum(lhs); sum += rhs; return sum;
}
template<unsigned nbits, unsigned rbits, typename bt>
constexpr tau<nbits, rbits, bt> operator-(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) {
	tau<nbits, rbits, bt> diff(lhs); diff -= rhs; return diff;
}
template<unsigned nbits, unsigned rbits, typename bt>
constexpr tau<nbits, rbits, bt> operator*(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) {
	tau<nbits, rbits, bt> mul(lhs); mul *= rhs; return mul;
}
template<unsigned nbits, unsigned rbits, typename bt>
constexpr tau<nbits, rbits, bt> operator/(const tau<nbits, rbits, bt>& lhs, const tau<nbits, rbits, bt>& rhs) {
	tau<nbits, rbits, bt> ratio(lhs); ratio /= rhs; return ratio;
}

// absolute value
template<unsigned nbits, unsigned rbits, typename bt>
constexpr tau<nbits, rbits, bt> abs(const tau<nbits, rbits, bt>& v) {
	return v.isneg() ? -v : v;
}

// components string for diagnostics
template<unsigned nbits, unsigned rbits, typename bt>
inline std::string components(const tau<nbits, rbits, bt>& v) {
	std::stringstream s;
	s << "(" << (v.sign() ? "-" : "+") << ", scale " << v.scale() << ", value " << v.to_double() << ")";
	return s.str();
}

}}  // namespace sw::universal

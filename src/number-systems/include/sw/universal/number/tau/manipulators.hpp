#pragma once
// manipulators.hpp: definitions of helper functions for balanced tau number manipulation
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <iostream>
#include <iomanip>
#include <sstream>
#include <typeinfo>  // for typeid()

#include <universal/native/integers.hpp>     // type_tag for native integer block types
#include <universal/utility/color_print.hpp> // color printing for shells

namespace sw { namespace universal {

	// Generate a type tag for this tau configuration
	template<typename TauType,
		std::enable_if_t< is_tau<TauType>, bool> = true
	>
	inline std::string type_tag(const TauType& = {}) {
		std::stringstream s;
		typename TauType::BlockType bt{ 0 };
		s << "tau<"
			<< std::setw(3) << TauType::nbits << ", "
			<< std::setw(3) << TauType::rbits << ", "
			<< type_tag(bt) << '>';
		return s.str();
	}

	// Generate a type field descriptor for this tau
	template<typename TauType,
		std::enable_if_t< is_tau<TauType>, bool> = true
	>
	inline std::string type_field(const TauType& = {}) {
		std::stringstream s;
		s << "fields(trits:" << TauType::nbits
			<< "|integer:" << TauType::ibits
			<< "|fraction:" << TauType::rbits << ')';
		return s.str();
	}

	// character used to print a balanced trit: '+' = +1, '0' = 0, '-' = -1
	inline char trit_char(int d) { return (d > 0) ? '+' : (d < 0 ? '-' : '0'); }

	// render the balanced-tau trit string, integer part '.' fraction part, prefix 0t
	template<unsigned nbits, unsigned rbits, typename bt>
	inline std::string to_triplets(const tau<nbits, rbits, bt>& v, bool tritMarker = false) {
		std::stringstream s;
		s << "0t";
		// integer trits: positions [nbits-1 .. rbits]
		if constexpr (nbits > rbits) {
			for (int p = static_cast<int>(nbits) - 1; p >= static_cast<int>(rbits); --p) {
				s << trit_char(v.trit(static_cast<unsigned>(p)));
				if (tritMarker && p > static_cast<int>(rbits) && ((p - static_cast<int>(rbits)) % 3) == 0) s << '\'';
			}
		}
		// fraction trits: positions [rbits-1 .. 0]
		if constexpr (rbits > 0) {
			s << '.';
			for (int p = static_cast<int>(rbits) - 1; p >= 0; --p) {
				s << trit_char(v.trit(static_cast<unsigned>(p)));
				if (tritMarker && p > 0 && (p % 3) == 0) s << '\'';
			}
		}
		return s.str();
	}

	// to_binary renders the underlying 2-bit-per-trit storage
	template<unsigned nbits, unsigned rbits, typename bt>
	inline std::string to_binary(const tau<nbits, rbits, bt>& v, bool nibbleMarker = false) {
		std::stringstream s;
		s << "0b";
		const auto& b = v.bits();
		for (int i = static_cast<int>(2 * nbits) - 1; i >= 0; --i) {
			s << (b.at(static_cast<unsigned>(i)) ? '1' : '0');
			if (nibbleMarker && i > 0 && (i % 4) == 0) s << '\'';
		}
		return s.str();
	}

	// native semantic representation: the balanced-tau trit string
	template<unsigned nbits, unsigned rbits, typename bt>
	inline std::string to_native(const tau<nbits, rbits, bt>& v, bool tritMarker = false) {
		return to_triplets(v, tritMarker);
	}

	template<unsigned nbits, unsigned rbits, typename bt>
	inline std::string components(const tau<nbits, rbits, bt>& v);  // defined in tau_impl.hpp

	template<typename TauType,
		std::enable_if_t< is_tau<TauType>, bool> = true
	>
	inline std::string range(const TauType& = {}) {
		std::stringstream s;
		TauType b(SpecificValue::maxneg), c(SpecificValue::minneg), d(SpecificValue::minpos), e(SpecificValue::maxpos);
		s << "[" << b << " ... " << c << ", 0, " << d << " ... " << e << "]";
		return s.str();
	}

	// report if a native floating-point value is within the dynamic range of the tau configuration
	template<typename TauType,
		std::enable_if_t< is_tau<TauType>, bool> = true
	>
	inline bool isInRange(double v) {
		TauType a{};
		return (v <= double(a.maxpos_bound()) && v >= -double(a.maxpos_bound()));
	}

	template<typename TauType,
		std::enable_if_t< is_tau<TauType>, bool> = true
	>
	inline std::string pretty_print(const TauType& v, bool tritMarker = false) {
		return to_triplets(v, tritMarker);
	}

	template<typename TauType,
		std::enable_if_t< is_tau<TauType>, bool> = true
	>
	inline std::string info_print(const TauType& = {}, int = 17) {
		return std::string("balanced base-tau (phi^2) number system");
	}

	template<unsigned nbits, unsigned rbits, typename bt>
	inline std::string color_print(const tau<nbits, rbits, bt>& v, bool = false) {
		std::stringstream s;
		Color cyan(ColorCode::FG_CYAN);
		Color magenta(ColorCode::FG_MAGENTA);
		Color def(ColorCode::FG_DEFAULT);
		if constexpr (nbits > rbits) {
			for (int p = static_cast<int>(nbits) - 1; p >= static_cast<int>(rbits); --p)
				s << cyan << trit_char(v.trit(static_cast<unsigned>(p)));
		}
		if constexpr (rbits > 0) {
			s << def << '.';
			for (int p = static_cast<int>(rbits) - 1; p >= 0; --p)
				s << magenta << trit_char(v.trit(static_cast<unsigned>(p)));
		}
		s << def;
		return s.str();
	}

}} // namespace sw::universal

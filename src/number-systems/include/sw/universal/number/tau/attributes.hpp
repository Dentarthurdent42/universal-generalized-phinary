#pragma once
// attributes.hpp: information functions for balanced tau number and value attributes
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <cmath> // for std::pow()

namespace sw { namespace universal {

// functions to provide details about the properties of a balanced tau number
// system configuration in terms of native types. Since many tau configurations
// cannot be represented by native types, these are convenience functions and
// should not be used for the core algorithms.

	// free function sign
template<unsigned nbits, unsigned rbits, typename bt>
bool sign(const tau<nbits, rbits, bt>& v) {
	return v.sign();
}

// free function scale
template<unsigned nbits, unsigned rbits, typename bt>
int scale(const tau<nbits, rbits, bt>& v) {
	return v.scale();
}

// generate the maxneg through maxpos value range of a tau configuration
template<unsigned nbits, unsigned rbits, typename bt>
std::string tau_range(const tau<nbits, rbits, bt>& v) {
	using TAU = tau<nbits, rbits, bt>;
	std::stringstream s;
	TAU l;
	s << std::setw(45) << type_tag(v) << " : [ "
		<< l.maxneg() << " ... "
		<< l.minneg() << " "
		<< "0 "
		<< l.minpos() << " ... "
		<< l.maxpos() << " ]";
	return s.str();
}

}} // namespace sw::universal

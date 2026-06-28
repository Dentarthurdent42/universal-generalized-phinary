#pragma once
// exceptions.hpp: definition of arbitrary configuration balanced tau number system exceptions
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <universal/common/exceptions.hpp>

namespace sw { namespace universal {

struct tau_arithmetic_exception : public universal_arithmetic_exception {
	tau_arithmetic_exception(const std::string& err) : universal_arithmetic_exception(std::string("tau arithmetic exception: ") + err) {}
};

// divide by zero arithmetic exception for reals
struct tau_divide_by_zero : public tau_arithmetic_exception {
	tau_divide_by_zero() : tau_arithmetic_exception("division by zero") {}
};

// negative argument to sqrt
struct tau_negative_sqrt_arg : public tau_arithmetic_exception {
	tau_negative_sqrt_arg() : tau_arithmetic_exception("negative sqrt argument") {}
};

///////////////////////////////////////////////////////////////
// internal implementation exceptions

// base class for internal exceptions
struct tau_internal_exception : public universal_internal_exception {
	tau_internal_exception(const std::string& error)
		: universal_internal_exception(std::string("tau internal error: ") + error) {};
};

struct tau_index_out_of_bounds : public tau_internal_exception {
	tau_index_out_of_bounds() : tau_internal_exception("index out of bounds") {}
};

}} // namespace sw::universal

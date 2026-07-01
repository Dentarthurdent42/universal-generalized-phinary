// tau.hpp: balanced tau (phi^2) number system standard header
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
//
// The balanced tau number system is the "ternary phi-representation" of the
// Generalized Phi Number System (GPNS) by S. Agaian and Y. Zhou,
// "Generalized Phi Number System and its Applications for Image Decomposition
// and Enhancement", Proc. SPIE-IS&T Electronic Imaging, Vol. 7881, 2011.
// It is a positional number system with the irrational base
//     tau = phi^2 = (3 + sqrt(5)) / 2
// and the balanced-ternary digit set { -1, 0, +1 }.
#pragma once
////////////////////////////////////////////////////////////////////////////////////////
///  COMPILATION DIRECTIVES TO DIFFERENT COMPILERS
#include <universal/utility/compiler.hpp>
#include <universal/utility/architecture.hpp>
#include <universal/utility/bit_cast.hpp>
#include <universal/utility/long_double.hpp>

////////////////////////////////////////////////////////////////////////////////////////
/// required std libraries
#include <iostream>
#include <iomanip>

////////////////////////////////////////////////////////////////////////////////////////
///  BEHAVIORAL COMPILATION SWITCHES

////////////////////////////////////////////////////////////////////////////////////////
// enable throwing specific exceptions for balanced tau arithmetic errors
// left to application to enable
#if !defined(TAU_THROW_ARITHMETIC_EXCEPTION)
// default is to use saturation/quiet handling instead of throwing
#define TAU_THROW_ARITHMETIC_EXCEPTION 0
#endif

///////////////////////////////////////////////////////////////////////////////////////
// bring in the trait functions
#include <universal/traits/number_traits.hpp>
#include <universal/traits/arithmetic_traits.hpp>

////////////////////////////////////////////////////////////////////////////////////////
/// INCLUDE FILES that make up the library
#include <universal/number/tau/exceptions.hpp>
#include <universal/number/tau/tau_fwd.hpp>
#include <universal/number/tau/tau_traits.hpp>
#include <universal/number/tau/tau_impl.hpp>
#include <universal/number/tau/numeric_limits.hpp>

// useful functions to work with tau numbers
#include <universal/number/tau/manipulators.hpp>
#include <universal/number/tau/attributes.hpp>

///////////////////////////////////////////////////////////////////////////////////////
/// math functions
#include <universal/number/tau/mathlib.hpp>

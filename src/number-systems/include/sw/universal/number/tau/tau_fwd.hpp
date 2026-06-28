#pragma once
// tau_fwd.hpp :  forward declarations of the balanced tau (phi^2) number system environment
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <cstddef>  // for unsigned

namespace sw { namespace universal {

// core tau types and functions
//   tau<nbits, rbits, bt> : balanced base-tau positional number system, tau = phi^2 = (3+sqrt(5))/2
//   nbits : total number of balanced-ternary trits (digits in {-1,0,+1})
//   rbits : number of fractional trits (positions -1 .. -rbits)
//   bt    : block storage type
template<unsigned nbits, unsigned rbits, typename bt> class tau;
template<unsigned nbits, unsigned rbits, typename bt> constexpr tau<nbits, rbits, bt> abs(const tau<nbits, rbits, bt>&);
template<unsigned nbits, unsigned rbits, typename bt> tau<nbits, rbits, bt> sqrt(const tau<nbits, rbits, bt>&);

}} // namespace sw::universal

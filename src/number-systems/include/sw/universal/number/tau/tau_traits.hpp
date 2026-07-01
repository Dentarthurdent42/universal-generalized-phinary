#pragma once
//  tau_traits.hpp : traits for the balanced tau (phi^2) number system
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
// SPDX-License-Identifier: MIT
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <universal/traits/integral_constant.hpp>

namespace sw { namespace universal {

// define a trait for tau types
template<typename _Ty>
struct is_tau_trait
	: false_type
{
};

template<unsigned nbits, unsigned rbits, typename BlockType>
struct is_tau_trait< tau<nbits, rbits, BlockType> >
	: true_type
{
};

template<typename _Ty>
constexpr bool is_tau = is_tau_trait<_Ty>::value;

template<typename _Ty>
using enable_if_tau = std::enable_if_t<is_tau<_Ty>, _Ty>;

}} // namespace sw::universal

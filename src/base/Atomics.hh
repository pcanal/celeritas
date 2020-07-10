//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file Atomics.hh
//---------------------------------------------------------------------------//
#ifndef base_Atomics_hh
#define base_Atomics_hh

#include "Assert.hh"
#include "Macros.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Add to a value, returning the original value.
 */
template<class T>
CELER_FORCEINLINE_FUNCTION T atomic_add(T* address, T value)
{
#ifdef __CUDA_ARCH__
    return atomicAdd(address, value);
#else
    REQUIRE(address);
    T initial = *address;
    *address += value;
    return initial;
#endif
}
//---------------------------------------------------------------------------//
} // namespace celeritas

#endif // base_Atomics_hh
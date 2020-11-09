//---------------------------------*-C++-*-----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SimStateStore.cc
//---------------------------------------------------------------------------//
#include "SimStateStore.hh"

#include "base/Array.hh"
#include "SimStatePointers.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with number of parallel tracks.
 */
SimStateStore::SimStateStore(size_type size) : vars_(size)
{
    REQUIRE(size > 0);
    ENSURE(!vars_.empty());
}

//---------------------------------------------------------------------------//
/*!
 * Number of threads stored in the state.
 */
size_type SimStateStore::size() const
{
    return vars_.size();
}

//---------------------------------------------------------------------------//
/*!
 * View to on-device state data.
 */
SimStatePointers SimStateStore::device_pointers()
{
    SimStatePointers result;
    result.vars = vars_.device_pointers();

    ENSURE(result);
    return result;
}

//---------------------------------------------------------------------------//
} // namespace celeritas

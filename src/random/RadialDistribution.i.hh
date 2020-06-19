//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file RadialDistribution.i.hh
//---------------------------------------------------------------------------//

#include <cmath>

#include "base/Assert.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with defaults.
 */
template<typename T>
CELER_INLINE_FUNCTION
RadialDistribution<T>::RadialDistribution(real_type radius)
    : radius_(radius), sample_uniform_{0, 1}
{
    REQUIRE(radius_ > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a random number according to the distribution
 */
template<class T>
template<class Generator>
CELER_INLINE_FUNCTION auto RadialDistribution<T>::operator()(Generator& rng)
    -> result_type
{
    return std::pow(this->sample_uniform_(rng),
                    static_cast<real_type>(1.0 / 3.0))
           * radius_;
}

//---------------------------------------------------------------------------//
} // namespace celeritas
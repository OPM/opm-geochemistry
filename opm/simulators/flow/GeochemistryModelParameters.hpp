/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef GEOCHEMISTRY_MODEL_PARAMETERS_HPP
#define GEOCHEMISTRY_MODEL_PARAMETERS_HPP

namespace Opm::Parameters
{

//! Target Courant number of the substeps of the explicit reactive transport
template <class Scalar>
struct GeochemistryTargetCfl
{
    static constexpr Scalar value = 1.0;
};

} // namespace Opm::Parameters

namespace Opm
{

/// Parameters of the geochemistry model
template <class Scalar>
struct GeochemistryModelParameters
{
    /// Read the runtime parameters, which must have been registered
    GeochemistryModelParameters();

    /// Register the runtime parameters
    static void registerParameters();

    /// Target Courant number of the substeps of the explicit reactive transport
    Scalar target_cfl_;
};

} // namespace Opm

#endif

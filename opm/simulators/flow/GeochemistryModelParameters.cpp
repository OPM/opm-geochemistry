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
#include <config.h>
#include <opm/simulators/flow/GeochemistryModelParameters.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>

namespace Opm
{

template <class Scalar>
GeochemistryModelParameters<Scalar>::GeochemistryModelParameters()
{
    target_cfl_ = Parameters::Get<Parameters::GeochemistryTargetCfl<Scalar>>();

    // The explicit transport is unstable for a Courant number above one, so the target is limited
    // to one. The number of substeps is also limited, so a very small target has little further
    // effect and is limited to a small positive value.
    constexpr Scalar minTargetCfl = 1.0e-2;
    constexpr Scalar maxTargetCfl = 1.0;
    const Scalar limited = std::isnan(target_cfl_)
        ? minTargetCfl
        : std::clamp(target_cfl_, minTargetCfl, maxTargetCfl);
    if (limited != target_cfl_) {
        OpmLog::warning(
            fmt::format("The target Courant number of the geochemistry transport is {}, "
                        "which is outside the range [{}, {}]. Using {}.",
                        target_cfl_,
                        minTargetCfl,
                        maxTargetCfl,
                        limited));
    }
    target_cfl_ = limited;
}

template <class Scalar>
void
GeochemistryModelParameters<Scalar>::registerParameters()
{
    Parameters::Register<Parameters::GeochemistryTargetCfl<Scalar>>(
        "Target Courant number for the explicit reactive transport. The number of "
        "reactive transport substeps is calculated to keep the Courant number below this value");
}

template struct GeochemistryModelParameters<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct GeochemistryModelParameters<float>;
#endif

} // namespace Opm

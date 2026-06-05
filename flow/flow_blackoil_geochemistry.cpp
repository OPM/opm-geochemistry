/*
  Copyright 2025, Equinor ASA

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

#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/FlowProblemGeochemistry.hpp>
#include <opm/simulators/flow/Main.hpp>

namespace Opm::Properties {
    namespace TTag {
        struct FlowProblemGeochem {
            using InheritsFrom = std::tuple<FlowProblem>;
        };
    }  // namespace TTag

    template<class TypeTag>
    struct Linearizer<TypeTag, TTag::FlowProblemGeochem> { using type = TpfaLinearizer<TypeTag>; };

    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::FlowProblemGeochem> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::FlowProblemGeochem> { static constexpr bool value = false; };

    template<class TypeTag>
    struct EnableGeochemistry<TypeTag, TTag::FlowProblemGeochem> { static constexpr bool value = true; };

    // Set the geochemistry problem
    template <class TypeTag>
    struct Problem<TypeTag, TTag::FlowProblemGeochem>
    {
        using type = FlowProblemGeochemistry<TypeTag>;
    };

}  // namspace Opm::Properties

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowProblemGeochem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}

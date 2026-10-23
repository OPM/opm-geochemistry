/*
  Copyright 2025 Equinor ASA.

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
#ifndef FLOW_PROBLEM_GEOCHEMISTRY_HPP
#define FLOW_PROBLEM_GEOCHEMISTRY_HPP

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/io/vtkgeochemistrymodule.hpp>
#include <opm/models/io/vtkgeochemistryparams.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/GeochemistryModel.hpp>

#include <memory>
#include <stdexcept>
#include <string>

namespace Opm {

template <typename TypeTag>
class FlowProblemGeochemistry : public FlowProblemBlackoil<TypeTag>
{
public:
    using Parent = FlowProblemBlackoil<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GeochemModel = GeochemistryModel<TypeTag>;

    enum { enableGeochem = getPropValue<TypeTag, Properties::EnableGeochemistry>() };

    /*!
    * \brief Constructor
    *
    * \param simulator Reference to simulator object
    */
    explicit FlowProblemGeochemistry(Simulator& simulator)
        : Parent(simulator)
        , geochemistryModel_(simulator)
    {
        const bool deckHasGeochem = simulator.vanguard().eclState().runspec().geochem().enabled();
        if constexpr (enableGeochem) {
            // Add VTK module
            this->model().addOutputModule(
                std::make_unique<VtkGeochemistryModule<TypeTag>>(simulator));

            // Sanity check
            if (!deckHasGeochem) {
                const std::string msg
                    = "Simulator with geochemistry enabled compile time, but deck "
                      "does not contain GEOCHEM keyword!";
                OpmLog::error(msg);
                throw std::runtime_error(msg);
            }
        } else {
            // Sanity check
            if (deckHasGeochem) {
                const std::string msg
                    = "GEOCHEM keyword in deck, but geochemistry disabled compile-time!";
                OpmLog::error(msg);
                throw std::runtime_error(msg);
            }
        }
    }

    /*!
    * \brief Register runtime parameters
    */
    static void registerParameters()
    {
        Parent::registerParameters();
        GeochemModel::registerParameters();
        VtkGeochemistryParams::registerParameters();
    }

    /*!
    * \brief Initialize the problem
    */
    void finishInit()
    {
        Parent::finishInit();
        geochemistryModel_.init();
    }

    /*!
    * \brief Calculations before a time integration
    */
    void beginTimeStep()
    {
        Parent::beginTimeStep();
        geochemistryModel_.beginTimeStep();
    }

    /*!
    * \brief Calculations after a time integration
    */
    void endTimeStep()
    {
        // Must come before Parent::endTimeStep(), which evaluates the summary state. Otherwise the
        // species rates in the summary lag one time step behind.
        geochemistryModel_.endTimeStep();
        Parent::endTimeStep();
    }

    /*!
    * \brief Get const reference to geochemistry model
    *
    * \returns Reference to geochemistry model
    */
    const GeochemModel& geochemistryModel() const
    {
        return geochemistryModel_;
    }

    /*!
    * \brief Get reference to geochemistry model
    *
    * \returns Reference to geochemistry model
    */
    GeochemModel& geochemistryModel()
    {
        return geochemistryModel_;
    }

    /*!
    * \brief Serialize variables
    *
    * \param serializer Byte array conversion
    */
    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<Parent&>(*this));
        serializer(geochemistryModel_);
    }

private:
    GeochemModel geochemistryModel_;
};  // class FlowProblemGeochemistry
}  // namespace Opm

#endif
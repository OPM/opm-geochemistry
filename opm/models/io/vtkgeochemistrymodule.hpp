// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef VTK_GEOCHEMISTRY_MODULE_HPP
#define VTK_GEOCHEMISTRY_MODULE_HPP

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkgeochemistryparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>


#include <string>
#include <vector>

namespace Opm {

template <class TypeTag>
class VtkGeochemistryModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    enum { enableGeochem = getPropValue<TypeTag, Properties::EnableGeochemistry>() };

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    /*!
    * \brief Constructor
    */
    explicit VtkGeochemistryModule(const Simulator& simulator)
            : ParentType(simulator)
    {
        if constexpr(enableGeochem) {
            params_.read();
        }
    }

    /*!
    * \brief Allocate output buffers
    */
    void allocBuffers() override
    {
        if constexpr(enableGeochem) {
            // Return if VTK is not enabled
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            // Initialize concentration output
            const auto& geochemistryModel = this->simulator_.problem().geochemistryModel();
            if (params_.speciesConcentrationOutput_) {
                speciesConcentration_.resize(geochemistryModel.numSpecies());
                for (std::size_t speciesIdx = 0; speciesIdx < geochemistryModel.numSpecies(); ++speciesIdx) {
                    this->resizeScalarBuffer_(speciesConcentration_[speciesIdx], BufferType::Dof);
                }
            }

            if (params_.mineralConcentrationOutput_ && geochemistryModel.numMinerals() > 0) {
                mineralConcentration_.resize(geochemistryModel.numMinerals());
                for (std::size_t mineralIdx = 0; mineralIdx < geochemistryModel.numMinerals(); ++mineralIdx) {
                    this->resizeScalarBuffer_(mineralConcentration_[mineralIdx], BufferType::Dof);
                }
            }
        }
    }

    /*!
    * \brief Process output variables
    *
    * \param elemCtx Reference to element context object
    */
    void processElement(const ElementContext& elemCtx) override
    {
        if constexpr(enableGeochem) {
            // Return if VTK is not enabled
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            // Get concentrations
            const auto& geochemistryModel = this->simulator_.problem().geochemistryModel();
            if (params_.speciesConcentrationOutput_) {
                for (std::size_t speciesIdx = 0; speciesIdx < geochemistryModel.numSpecies(); ++speciesIdx) {
                    for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                        const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                        speciesConcentration_[speciesIdx][globalDofIdx] =
                            geochemistryModel.speciesConcentration(speciesIdx, globalDofIdx);
                    }
                }
            }

            if (params_.mineralConcentrationOutput_ && geochemistryModel.numMinerals() > 0) {
                for (std::size_t mineralIdx = 0; mineralIdx < geochemistryModel.numMinerals(); ++mineralIdx) {
                    for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                        const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                        mineralConcentration_[mineralIdx][globalDofIdx] =
                            geochemistryModel.mineralConcentration(mineralIdx, globalDofIdx);
                    }
                }
            }
        }
    }

    /*!
    * \brief Add output to file
    *
    * \param baseWriter Reference to output writer
    */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if constexpr(enableGeochem) {
            // Return if VTK is not enabled
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            // Return if basewriter does not exist
            if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
                return;
            }

            // Commit concentration output
            const auto& geochemistryModel = this->simulator_.problem().geochemistryModel();
            if (params_.speciesConcentrationOutput_) {
                for (std::size_t speciesIdx = 0; speciesIdx < geochemistryModel.numSpecies(); ++speciesIdx) {
                    const std::string cname = "speciesConcentration_" + geochemistryModel.speciesName(speciesIdx);
                    this->commitScalarBuffer_(baseWriter,
                                              cname.c_str(),
                                              speciesConcentration_[speciesIdx],
                                              BufferType::Dof);
                }

            }
            if (params_.mineralConcentrationOutput_ && geochemistryModel.numMinerals() > 0) {
                for (std::size_t mineralIdx = 0; mineralIdx < geochemistryModel.numMinerals(); ++mineralIdx) {
                    const std::string mname = "mineralConcentration_" + geochemistryModel.mineralName(mineralIdx);
                    this->commitScalarBuffer_(baseWriter,
                                              mname.c_str(),
                                              mineralConcentration_[mineralIdx],
                                              BufferType::Dof);
                }
            }
        }
    }

private:
    VtkGeochemistryParams params_{};
    std::vector<ScalarBuffer> speciesConcentration_;
    std::vector<ScalarBuffer> mineralConcentration_;
};  // class VtkGeochemistryModule

} // namespace Opm

#endif
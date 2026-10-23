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
#ifndef GEOCHEMISTRY_MODEL_HPP
#define GEOCHEMISTRY_MODEL_HPP

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/bvector.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Geochemistry/GenericSpeciesConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>

#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/GeochemistryModelParameters.hpp>
#include <opm/simulators/geochemistry/OpmGeoChemInterface.hpp>
#include <opm/simulators/wells/WellTracerRate.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm {

template <class TypeTag>
class GeochemistryModel
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using SpeciesVector = Dune::BlockVector<Scalar>;

    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    static constexpr int maxNumSubsteps = 100;
    static constexpr Scalar tolSat = 1e-6; //!< Water saturation below which a cell is dry
    static constexpr Scalar minSatInitial = 0.1; //!< Minimum water saturation of initial solve

    struct WellSource
    {
        std::vector<Scalar> speciesInjectionRate; //!< Species amount injected per time
        Scalar waterProductionRate {0}; //!< Water produced per time
    };
    struct ProducerConnection
    {
        std::size_t wellSeqIndex; //!< Sequence index of the well
        std::size_t cell; //!< Cell of the connection
        Scalar rate; //!< Water rate of the connection
        int segment; //!< Segment of the connection for a multisegment well, otherwise -1
    };
    struct FaceFlux
    {
        unsigned neighbour; //!< Cell on the other side of the face
        Scalar flux; //!< Water flux out of the cell, in surface volume per time
        bool isUp; //!< True if the cell itself is upstream
    };

public:
    /*!
    * \brief Constructor
    *
    * \param simulator Reference to simulator object
    */
    explicit GeochemistryModel(Simulator& simulator)
        : simulator_ (simulator)
        , eclState_  (simulator.vanguard().eclState())
        , cartMapper_(simulator.vanguard().cartesianIndexMapper())
        , element_chunks_(simulator.gridView(), Dune::Partitions::all, ThreadManager::maxThreads())
    {}

    /*!
     * \brief Register runtime parameters
     */
    static void registerParameters()
    {
        GeochemistryModelParameters<Scalar>::registerParameters();
    }

    /*!
    * \brief Initialize geochemistry model
    *
    * Initialize the geochemistry solver and concentration vectors
    */
    void init()
    {
        // Return if GEOCHEM is not activated
        const auto& geochem = eclState_.runspec().geochem();
        if (!geochem.enabled()) {
            return;
        }

        // Get transported species names
        const auto& species = eclState_.species();
        speciesNames_.reserve(species.size());
        std::ranges::transform(species, std::back_inserter(speciesNames_),
                               [] (const auto& item) { return item.name; } );

        // Get mineral names
        const auto& mineral = eclState_.mineral();
        if (!mineral.empty()) {
            mineralNames_.reserve(mineral.size());
            std::ranges::transform(mineral, std::back_inserter(mineralNames_),
                                   [] (const auto& item) { return item.name; } );
        }

        // Get ion exchange names
        const auto& ion_exchange = eclState_.ionExchange();
        if (!ion_exchange.empty()) {
            ionExNames_.reserve(ion_exchange.size());
            std::ranges::transform(ion_exchange, std::back_inserter(ionExNames_),
                                   [] (const auto& item) { return item.name; } );
        }

        // Initialize interface to geochemistry solver
        const auto& file_name = geochem.geochem_file_name();
        std::pair<double, double> tol =
                std::make_pair<double, double>(geochem.mbal_tol(), geochem.ph_tol());
        bool charge_balance = geochem.charge_balance();
        int splay_tree_resolution = geochem.splay_tree_resolution();
        geoChemInterface_ = std::make_shared<OpmGeoChemInterface>();
        geoChemInterface_->initialize_from_opm_deck(file_name,
                                                    speciesNames_,
                                                    mineralNames_,
                                                    ionExNames_,
                                                    charge_balance,
                                                    tol,
                                                    splay_tree_resolution);

        // Minerals and ion exchangers must be the same in the deck and the geochemistry solver
        if (geoChemInterface_->numberOfMinerals() != numMinerals()) {
            throw std::runtime_error(
                fmt::format("The geochemistry solver has {} minerals, but {} are given by the "
                            "MINERAL keyword in the deck.",
                            geoChemInterface_->numberOfMinerals(),
                            numMinerals()));
        }
        // The solver orders minerals differently from the deck, hence an index mapping is needed
        deckToSolverMineral_ = geoChemInterface_->mineral_index_map(mineralNames_);
        if (geoChemInterface_->numberOfIonExchange() != numIonEx()) {
            throw std::runtime_error(
                fmt::format("The geochemistry solver has {} ion exchange species, but {} are given "
                            "by the IONEX keyword in the deck.",
                            geoChemInterface_->numberOfIonExchange(),
                            numIonEx()));
        }

        // Initialize species independent vectors
        const std::size_t numGridDof = simulator_.model().numGridDof();
        pH_.resize(numGridDof, 7.0);
        sigma_.resize(numGridDof, 0.0);
        psi_.resize(numGridDof, 0.0);
        vol1_.resize(numGridDof);
        volumeNew_.resize(numGridDof);
        faceFluxes_.resize(numGridDof);

        // Fill in species concentrations
        const std::size_t nSpecies = numSpecies();
        concentration_.resize(nSpecies);
        concentrationInitial_.resize(nSpecies);
        Cads_.resize(nSpecies);
        concentrationNext_.resize(nSpecies);
        for (std::size_t speciesIdx = 0; speciesIdx < nSpecies;  ++speciesIdx) {
            const auto& single_species = species[speciesIdx];
            concentration_[speciesIdx].resize(numGridDof);
            concentrationInitial_[speciesIdx].resize(numGridDof);
            Cads_[speciesIdx].resize(numGridDof);
            concentrationNext_[speciesIdx].resize(numGridDof);

            // Initial concentration for species
            setInitialConcentrations_(single_species, concentration_[speciesIdx]);
        }

        // Initialize mineral concentrations
        const std::size_t nMin = geoChemInterface_->numberOfMinerals();
        if (nMin > 0) {
            Cmin_.resize(nMin);
            minWt_.resize(nMin);
            for (std::size_t minSpeciesIdx = 0; minSpeciesIdx < nMin; ++minSpeciesIdx) {
                const auto& single_mineral = mineral[minSpeciesIdx];
                Cmin_[minSpeciesIdx].resize(numGridDof);
                minWt_[minSpeciesIdx].resize(numGridDof);

                // Initial weight fraction for mineral
                setInitialConcentrations_(single_mineral, minWt_[minSpeciesIdx]);
            }
        }

        // Initialize ion exchange species
        const std::size_t nIo = geoChemInterface_->numberOfIonExchange();
        if (nIo > 0) {
            Cio_.resize(nIo);
            for (std::size_t ioIdx = 0; ioIdx < nIo; ++ioIdx) {
                const auto& single_io = ion_exchange[ioIdx];
                Cio_[ioIdx].resize(numGridDof);

                // Initial ion exchange
                setInitialConcentrations_(single_io, Cio_[ioIdx]);
            }
        }
    }

    /*!
    * \brief Calculations before a time integration
    *
    * Updates variables from the previous time step to be used in endTimeStep()
    */
    void beginTimeStep()
    {
        // Return if GEOCHEM is not activated
        if (!eclState_.runspec().geochem().enabled()) {
            return;
        }

        // Store variables from previous time step
        updateStorageCache();

        // Equilibrate the initial state, which is only done in the first time step
        if (!initialEquilibrationDone_) {
            initialEquilibration_();
        }
    }

    /*!
    * \brief Calculations after a time integration
    *
    * \note The reactive transport step is done here!
    */
    void endTimeStep()
    {
        if (!eclState_.runspec().geochem().enabled()) {
            return;
        }

        // Reactive transport
        advanceSpeciesFieldsExplicit();
    }

    /*!
    * \brief Get number of species
    *
    * \returns Number of transported species
    */
    std::size_t numSpecies() const
    {
        return speciesNames_.size();
    }

    /*!
    * \brief Get number of minerals
    *
    * \returns Number of minerals
    */
    std::size_t numMinerals() const
    {
        return mineralNames_.size();
    }

    /*!
    * \brief Get number of ion exchange species
    *
    * \returns Number of ion exchange
    */
    std::size_t numIonEx() const
    {
        return ionExNames_.size();
    }

    /*!
    * \brief Get a particular species name
    *
    * \param speciesIdx Index of species
    * \returns String with queried species name
    */
    const std::string& speciesName(unsigned speciesIdx) const
    {
        return speciesNames_[speciesIdx];
    }

    /*!
    * \brief Get (reference to) a particular mineral name
    *
    * \param minSpeciesIdx Index of mineral
    * \returns String with queried mineral name
    */
    const std::string& mineralName(unsigned minSpeciesIdx) const
    {

        return mineralNames_[minSpeciesIdx];
    }

    /*!
    * \brief Get (reference to) a particular ion exchange name
    *
    * \param minSpeciesIdx Index of ion exchange
    * \returns String with queried ion exchage name
    */
    const std::string& ionExchangeName(unsigned ionExIdx) const
    {

        return ionExNames_[ionExIdx];
    }

    /*!
    * \brief Get concentration of a species at a grid cell
    *
    * \param speciesIdx Index of species
    * \param globalDofIdx Cell index
    * \returns Concentration of species
    */
    Scalar speciesConcentration(int speciesIdx, int globalDofIdx) const
    {
        if (concentration_.empty()) {
            return 0.0;
        }

        return concentration_[speciesIdx][globalDofIdx];
    }

    /*!
    * \brief Get concentration of a mineral at a grid cell
    *
    * \param minSpeciesIdx Index of mineral (deck order)
    * \param globalDofIdx Cell index
    * \returns Concentration of mineral
    */
    Scalar mineralConcentration(int minSpeciesIdx, int globalDofIdx) const
    {
        if (Cmin_.empty()) {
            return 0.0;
        }

        // Cmin_ is in solver order, minSpeciesIdx is in deck order
        return Cmin_[deckToSolverMineral_[minSpeciesIdx]][globalDofIdx];
    }

    /*!
    * \brief Get pH in a grid cell
    *
    * \param globalDofIdx Cell index
    * \returns pH
    */
    Scalar PH(int globalDofIdx) const
    {
        if (pH_.empty()) {
            return 0.0;
        }

        return pH_[globalDofIdx];
    }

    /*!
    * \brief Get all "standard" wells' species concentration rates
    *
    * \returns Container with well species concentration rates
    */
    const std::unordered_map<int, std::vector<WellTracerRate<Scalar>>>&
    getWellSpeciesRates() const
    {
        return wellSpeciesRate_;
    }

    /*!
    * \brief Get all multisegmented wells' species concentration rates
    *
    * \returns Container with well species concentration rates
    */
    const std::unordered_map<int, std::vector<MSWellTracerRate<Scalar>>>&
    getMswSpeciesRates() const
    {
        return mSwSpeciesRate_;
    }

    /*!
    * \brief Set species contration rate in a grid cell
    *
    * \param speciesIdx Index of mineral
    * \param globalDofIdx Cell index
    * \param value Species concentration for grid cell
    */
    void setSpeciesConcentration(int speciesIdx, int globalDofIdx, Scalar value)
    {
        concentration_[speciesIdx][globalDofIdx] = value;
    }

    /*!
    * \brief Serialize variables
    *
    * \param serializer Byte array conversion
    */
    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(concentration_);
        serializer(Cads_);
        serializer(Cmin_);
        serializer(pH_);
        serializer(sigma_);
        serializer(psi_);
        serializer(initialEquilibrationDone_);
        serializer(wellSpeciesRate_);
        serializer(mSwSpeciesRate_);
    }


protected:
    /*!
     * \brief Set initial concentrations for a (single) species
     *
     * \param single_species Species config
     * \param concentration Initial concentration vector to be set
     */
    void setInitialConcentrations_(const GenericSpeciesConfig::SpeciesEntry& single_species,
                                   SpeciesVector& concentration)
    {
        // *BLK
        if (single_species.concentration.has_value()) {
            const auto& species_concentration = single_species.concentration.value();
            if (species_concentration.size() != static_cast<std::size_t>(cartMapper_.cartesianSize())) {
                throw std::runtime_error("Size of S/M/IBLK " + single_species.name + " is wrong!");
            }

            // The *BLK concentrations are given for each Cartesian cell, while the concentration
            // vector has an entry for each active cell
            for (std::size_t globalDofIdx = 0; globalDofIdx < concentration.size();
                 ++globalDofIdx) {
                const int cartDofIdx = cartMapper_.cartesianIndex(globalDofIdx);
                concentration[globalDofIdx] = species_concentration[cartDofIdx];
            }
        }
        // *VDP
        else if (single_species.svdp.has_value()) {
            const auto& species_svdp = single_species.svdp.value();
            const auto& centroids = simulator_.vanguard().cellCentroids();

            // For each each grid cell, evaluate *VDP and assign to concentration vector
            std::for_each(
                concentration.begin(), concentration.end(),
                [idx = 0, &species_svdp, &centroids](auto& conc) mutable
                {
                    conc = species_svdp.evaluate("SPECIES_CONCENTRATION", centroids(idx)[2]);
                    ++idx;
                }
            );
        }
        // Zero initial condition
        else {
            OpmLog::warning(fmt::format("No S/M/IBLK or S/M/IVDP given for species {}. "
                                        "Initial values set to zero. ", single_species.name));
            std::ranges::fill(concentration, 0.0);
        }
    }

    /*!
    * \brief Compute black oil equation volume term
    *
    * \param globalDofIdx Cell index
    * \param timeIdx Time step index
    * \return Black oil volume term
    */
    Scalar computeVolume_(const unsigned globalDofIdx,
                          const unsigned timeIdx) const
    {
        const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        constexpr Scalar min_volume = 1e-10;

        return std::max(decay<Scalar>(fs.saturation(waterPhaseIdx)) *
                        decay<Scalar>(fs.invB(waterPhaseIdx)) *
                        decay<Scalar>(intQuants.porosity()),
                        min_volume);
    }

    /*!
    * \brief Compute black oil equation flux term
    *
    * \param elemCtx Reference to element context object
    * \param scvfIdx (Local) control volume index
    * \param timeIdx Time step index
    * \return Black oil flux term and boolean indicating upstream cell or not
    */
    std::pair<Scalar, bool> computeFlux_(const ElementContext& elemCtx,
                                         const unsigned scvfIdx,
                                         const unsigned timeIdx) const
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const unsigned inIdx = extQuants.interiorIndex();

        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        Scalar v = decay<Scalar>(extQuants.volumeFlux(waterPhaseIdx))
                   * decay<Scalar>(fs.invB(waterPhaseIdx));

        const Scalar A = scvf.area();
        return std::pair{A * v, inIdx == upIdx};
    }

    /*!
    * \brief Update cache for reactive transport step
    *
    * \warning Everything that is need from Flow -and is not stored before endTimeStep()- must be saved here!
    */
    void updateStorageCache()
    {
        // Update concentration from previous time step
        concentrationInitial_ = concentration_;

        // Parallel loop over element chunks
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (const auto& chunk : element_chunks_) {
            ElementContext elemCtx(simulator_);

            for (const auto& elem : chunk) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const Scalar extrusionFactor = elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
                const Scalar scvVolume = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume() * extrusionFactor;
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(0, /*timeIdx=*/0);

                vol1_[globalDofIdx] = computeVolume_(globalDofIdx, 0) * scvVolume;
            }
        }
    }

    /*!
     * \brief Send the concentrations in the cells owned by this rank to the ranks that have them as
     *        overlap cells
     *
     * \param concentrations Concentration vector for each species
     *
     * \note VectorVectorDataHandle.hpp has no include guard, so it is not included here again. It
     * is included by TracerModel.hpp, which FlowProblem.hpp includes before this file.
     */
    void communicateConcentrations_(std::vector<SpeciesVector>& concentrations)
    {
        auto handle = VectorVectorDataHandle<GridView, std::vector<SpeciesVector>>(
            concentrations, simulator_.gridView());
        simulator_.gridView().communicate(
            handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
    }

    /*!
    * \brief Equilibrate the geochemical system in all interior cells with the initial concentrations
    *
    * \note Done once, in the first beginTimeStep(), so every cell is equilibrated before any
    * transport is calculated. Uses the concentrations stored by updateStorageCache().
    */
    void initialEquilibration_()
    {
        ElementContext elemCtx(simulator_);
        for (const auto& elem : elements(simulator_.gridView())) {
            if (elem.partitionType() != Dune::InteriorEntity) {
                continue;
            }

            elemCtx.updateStencil(elem);
            const unsigned I = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/0);
            speciesEquationChemistryExplicit_(/*dt=*/0.0, I, /*initial_equil=*/true);
        }

        // The overlap cells, which are owned by other ranks, must have the equilibrated values too
        communicateConcentrations_(concentrationInitial_);

        initialEquilibrationDone_ = true;
    }

    /*!
    * \brief Run reactive transport solver and post-processing
    *
    * \note The actual reactive transport solver is in speciesEquationsExplicit_()
    */
    void advanceSpeciesFieldsExplicit()
    {
        // Calculate new concentration fields for each species
        speciesEquationsExplicit_();

        // Post-processing for concentration output
        // Only interior cells are processed, the overlap cells are overwritten by the
        // communication below.
        for (const auto globalDofIdx : interiorDofs_) {
            const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, 0);
            const auto& fs = intQuants.fluidState();
            const Scalar Sw = decay<Scalar>(fs.saturation(waterPhaseIdx));

            for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                if (concentration_[sIdx][globalDofIdx] < 0.0 || Sw < tolSat) {
                    concentration_[sIdx][globalDofIdx] = 0.0;
                }
            }
        }

        // The overlap cells are not calculated here, they get the values from the rank that owns
        // them. This is needed for the flux terms of the next time step.
        communicateConcentrations_(concentration_);

        // Correct the reported rates of the producers with cross flow, then convert to raw rates
        correctCrossFlowRates_();
        convertEffectiveRatesToRawRates_();
    }

    /*!
     * \brief Correct the reported species rates of the producers with cross flow
     *
     * If the well rate is larger than the sum of the producing connection rates, some connections
     * inject. The reported rates are then scaled with the ratio of the well rate to the sum of
     * the producing rates, or set to zero if the well rate is below a small threshold.
     */
    void correctCrossFlowRates_()
    {
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        for (const auto& wellPtr : wellPtrs) {
            const auto& eclWell = wellPtr->wellEcl();
            if (!eclWell.isProducer()) {
                continue;
            }

            const std::size_t well_index
                = simulator_.problem().wellModel().wellState().index(eclWell.name()).value();
            const auto& ws = simulator_.problem().wellModel().wellState().well(well_index);

            Scalar rateWellNeg = 0.0;
            for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
                const auto I = ws.perf_data.cell_index[i];
                const Scalar rate = wellPtr->volumetricSurfaceRateForConnection(I, waterPhaseIdx);
                if (rate < 0.0) {
                    rateWellNeg += rate;
                }
            }

            // TODO: Some inconsistencies here that perhaps should be clarified. The "official" rate
            // is occasionally significantly different from the sum over connections. Only observed
            // for small values, negligible for the rate itself, but matters when used to calculate
            // species concentrations.
            const Scalar rateWellTotal = ws.surface_rates[waterPhaseIdx];

            // The well rate and the connection rates differ by the tolerance of the well solver
            // also without cross flow, so a small relative difference is not cross flow.
            constexpr Scalar crossFlowTolerance = 1.0e-6;
            if (rateWellTotal - rateWellNeg
                > crossFlowTolerance * std::abs(rateWellNeg)) { // Cross flow
                constexpr Scalar bucketPrDay
                    = 10.0 / (1000. * 3600. * 24.); // ... keeps (some) trouble away
                const Scalar factor
                    = (rateWellTotal < -bucketPrDay) ? rateWellTotal / rateWellNeg : 0.0;
                for (auto& speciesRate : wellSpeciesRate_[eclWell.seqIndex()]) {
                    speciesRate.rate *= factor;
                }
            }
        }
    }

    /*!
     * \brief Convert the reported well species rates from effective to raw rates
     *
     * The connection rates from the well model include the well efficiency factor, which is right
     * for the species that are injected and produced in the cells. The reported rates should not
     * include it.
     */
    void convertEffectiveRatesToRawRates_()
    {
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        for (const auto& wellPtr : wellPtrs) {
            const auto& eclWell = wellPtr->wellEcl();
            const auto wellSeqIndex = eclWell.seqIndex();
            const auto invWellEffFactor
                = 1.0 / std::max<Scalar>(1.0e-10, wellPtr->wellEfficiencyFactor());

            std::ranges::for_each(wellSpeciesRate_[wellSeqIndex], [&](WellTracerRate<Scalar>& wsr) {
                wsr.rate *= invWellEffFactor;
            });
            if (eclWell.isMultiSegment()) {
                std::ranges::for_each(
                    mSwSpeciesRate_[wellSeqIndex], [&](MSWellTracerRate<Scalar>& wsr) {
                        std::ranges::for_each(wsr.rate,
                                              [&](auto& item) { item.second *= invWellEffFactor; });
                    });
            }
        }
    }

    /*!
     * \brief Reactive transport solver with explicit scheme
     *
     * The time step is split in as many substeps as needed to keep the Courant number below 1, so
     * that no cell loses more water than it holds in a substep. The fluxes and well rates are those
     * of the time step, and the pore volume changes linearly from the old to the new value. Each
     * substep does the transport and then the chemistry, so the chemistry follows the water as it
     * moves through the cells. With a Courant number below 1 there is a single substep.
     */
    void speciesEquationsExplicit_()
    {
        // Clear well containers
        wellSpeciesRate_.clear();
        mSwSpeciesRate_.clear();

        // Reserve new space
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        wellSpeciesRate_.reserve(wellPtrs.size());
        mSwSpeciesRate_.reserve(wellPtrs.size());

        // Simulator information
        ElementContext elemCtx(simulator_);
        const Scalar dt = elemCtx.simulator().timeStepSize();

        // Calculate well terms. The perforations of a rank are only in its interior cells.
        wellSources_.clear();
        producerConnections_.clear();
        for (const auto& wellPtr : wellPtrs) {
            speciesEquationWellExplicit_(*wellPtr);
        }

        // Quantities that are constant during the time step
        prepareTransportExplicit_(elemCtx);

        // Transport from the concentrations at the start of the time step
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            concentration_[sIdx] = concentrationInitial_[sIdx];
        }
        const int numSubsteps = numberOfSubsteps_(dt);
        for (int substep = 0; substep < numSubsteps; ++substep) {
            transportSubstepExplicit_(substep, numSubsteps, dt / numSubsteps);

            // Equilibrate geochemical system
            const Scalar substepFraction = static_cast<Scalar>(substep + 1) / numSubsteps;
            for (const auto I : interiorDofs_) {
                speciesEquationChemistryExplicit_(dt / numSubsteps, I, /*initial_equil=*/false,
                                                  substepFraction);
            }
            communicateConcentrations_(concentration_);
        }
    }

    /*!
     * \brief Calculate the quantities of the interior cells that are constant during a time step
     *
     * \param elemCtx Reference to element context object
     */
    void prepareTransportExplicit_(ElementContext& elemCtx)
    {
        interiorDofs_.clear();
        for (const auto& elem : elements(simulator_.gridView())) {
            if (elem.partitionType() != Dune::InteriorEntity) {
                continue;
            }

            elemCtx.updateStencil(elem);
            const unsigned I = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
            interiorDofs_.push_back(I);

            // Update block quantities
            elemCtx.updateAllIntensiveQuantities();
            elemCtx.updateAllExtensiveQuantities();

            // Volume at current time step
            const Scalar extrusionFactor
                = elemCtx.intensiveQuantities(/*dofIdx=*/0, /*timeIdx=*/0).extrusionFactor();
            Valgrind::CheckDefined(extrusionFactor);
            assert(isfinite(extrusionFactor));
            assert(extrusionFactor > 0.0);
            const Scalar scvVolume
                = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/0).volume()
                * extrusionFactor;
            volumeNew_[I] = computeVolume_(I, 0) * scvVolume;

            // Fluxes over the faces to the neighbour cells
            faceFluxes_[I].clear();
            const std::size_t numInteriorFaces = elemCtx.numInteriorFaces(/*timIdx=*/0);
            for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
                const auto& face = elemCtx.stencil(0).interiorFace(scvfIdx);
                const unsigned j = face.exteriorIndex();
                const unsigned J = elemCtx.globalSpaceIndex(/*dofIdx=*/j, /*timIdx=*/0);
                const auto& [flux, isUp] = computeFlux_(elemCtx, scvfIdx, 0);
                faceFluxes_[I].push_back({J, flux, isUp});
            }
        }
    }

    /*!
     * \brief Number of transport substeps needed to keep the Courant number below the target
     *
     * The Courant number of a cell is the water that leaves the cell in the time step (over its
     * faces and through wells) divided by the water the cell holds. The substeps are as few as
     * needed to keep the Courant number of each cell below the target Courant number. All ranks
     * must use the same number, as the overlap cells are communicated in each substep.
     *
     * \param dt Time step
     * \returns Number of substeps
     */
    int numberOfSubsteps_(const Scalar dt) const
    {
        // Water leaving the cell over its faces and through its wells, divided by its water volume
        // and by the target Courant number, so a value of one means that the cell is at the target
        auto courantRatio = [this, dt](const unsigned I, const Scalar wellOutflow) {
            Scalar outflow = wellOutflow;
            for (const auto& face : faceFluxes_[I]) {
                outflow += std::max<Scalar>(face.flux, 0);
            }
            return dt * outflow / (param_.target_cfl_ * std::min(vol1_[I], volumeNew_[I]));
        };

        Scalar maxCourantRatio = 0.0;
        for (const auto I : interiorDofs_) {
            maxCourantRatio = std::max(maxCourantRatio, courantRatio(I, 0));
        }
        for (const auto& [I, source] : wellSources_) {
            maxCourantRatio = std::max(
                maxCourantRatio, courantRatio(I, std::max<Scalar>(-source.waterProductionRate, 0)));
        }

        // All ranks must use the same number of substeps
        maxCourantRatio = simulator_.gridView().comm().max(maxCourantRatio);

        // Compare as Scalar, as converting a value above INT_MAX to int is undefined
        const Scalar neededSubsteps = std::ceil(maxCourantRatio);
        int numSubsteps = 1;
        if (neededSubsteps > static_cast<Scalar>(maxNumSubsteps)) {
            OpmLog::warning(
                fmt::format("The geochemistry transport needs {:.3g} substeps in a time step "
                            "of {} s, but is limited to {}. Concentrations can become "
                            "negative and are then set to zero. Reduce the time step.",
                            neededSubsteps,
                            dt,
                            maxNumSubsteps));
            numSubsteps = maxNumSubsteps;
        } else {
            numSubsteps = std::max(1, static_cast<int>(neededSubsteps));
            if (numSubsteps > 1) {
                OpmLog::debug(
                    fmt::format("Geochemistry transport uses {} substeps in a time step of {} s",
                                numSubsteps,
                                dt));
            }
        }
        return numSubsteps;
    }

    /*!
     * \brief One explicit upwind substep of the transport of the species
     *
     * \param substep Index of the substep
     * \param numSubsteps Number of substeps
     * \param dtSubstep Length of the substep
     */
    void transportSubstepExplicit_(const int substep, const int numSubsteps, const Scalar dtSubstep)
    {
        // Volume and flux terms
        const Scalar fractionStart = static_cast<Scalar>(substep) / numSubsteps;
        const Scalar fractionEnd = static_cast<Scalar>(substep + 1) / numSubsteps;
        for (const auto I : interiorDofs_) {
            // Volumes at start and end of substep (linear fraction of total volume change)
            const Scalar volumeStart = vol1_[I] + fractionStart * (volumeNew_[I] - vol1_[I]);
            const Scalar volumeEnd = vol1_[I] + fractionEnd * (volumeNew_[I] - vol1_[I]);
            for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                // Amount at the start and the fluxes from the concentrations at the start
                Scalar amount = volumeStart * concentration_[sIdx][I];
                for (const auto& face : faceFluxes_[I]) {
                    const unsigned upstream = face.isUp ? I : face.neighbour;
                    amount -= dtSubstep * face.flux * concentration_[sIdx][upstream];
                }
                concentrationNext_[sIdx][I] = amount / volumeEnd;
            }
        }

        // Injection/production from wells (loop only over well cells)
        for (const auto& [I, source] : wellSources_) {
            const Scalar volumeEnd = vol1_[I] + fractionEnd * (volumeNew_[I] - vol1_[I]);
            for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                const Scalar rate = source.speciesInjectionRate[sIdx]
                    + source.waterProductionRate * concentration_[sIdx][I];
                concentrationNext_[sIdx][I] += dtSubstep * rate / volumeEnd;
            }
        }

        // Report the species produced by the producers, which is the amount removed above, averaged
        // over the time step. The concentrations are still those at the start of the substep.
        const Scalar substepWeight = 1.0 / numSubsteps;
        for (const auto& connection : producerConnections_) {
            auto& speciesRate = wellSpeciesRate_[connection.wellSeqIndex];
            for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                const Scalar delta
                    = substepWeight * connection.rate * concentration_[sIdx][connection.cell];
                speciesRate[sIdx].rate += delta;
                if (connection.segment >= 0) {
                    mSwSpeciesRate_[connection.wellSeqIndex][sIdx].rate[connection.segment]
                        += delta;
                }
            }
        }

        // All cells are updated from the concentrations at the start, then the new values are used
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            for (const auto I : interiorDofs_) {
                concentration_[sIdx][I] = concentrationNext_[sIdx][I];
            }
        }
    }

    /*!
    * \brief Get WSPECIES for a particular species
    *
    * \param eclWell Reference to eclWell object
    * \param name Name of species
    * \param summaryState SummaryState object
    * \returns Concentration of injected species
    */
    Scalar currentWSPECIES_(const Well& eclWell, const std::string& name, const SummaryState& summaryState) const
    {
        return eclWell.getSpeciesProperties().getConcentration(WellTracerProperties::Well { eclWell.name() },
                                                               WellTracerProperties::Tracer { name },
                                                               summaryState);
    }

    /*!
     * \brief Get the well source of a cell, which is created if the cell has none yet
     *
     * \param I Cell index
     */
    WellSource& getWellSource_(const unsigned I)
    {
        auto& source = wellSources_[I];
        if (source.speciesInjectionRate.empty()) {
            source.speciesInjectionRate.assign(numSpecies(), 0);
        }
        return source;
    }

    /*!
     * \brief Calculate single well contribution in explicit reactive transport solver
     *
     * The rates are per time. They are used in each transport substep.
     *
     * \param well Reference to well object
     */
    template <class Well>
    void speciesEquationWellExplicit_(const Well& well)
    {
        // Get simulation wells
        const auto& eclWell = well.wellEcl();

        // Reserve space for species output
        auto& speciesRate = wellSpeciesRate_[eclWell.seqIndex()];
        speciesRate.reserve(numSpecies());
        const bool isMsw = eclWell.isMultiSegment();
        auto* mswSpeciesRate = isMsw
            ? &mSwSpeciesRate_[eclWell.seqIndex()]
            : nullptr;
        if (mswSpeciesRate) {
            mswSpeciesRate->reserve(numSpecies());
        }

        // Init. well output to zero
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            speciesRate.emplace_back(speciesName(sIdx), 0.0);
            if (isMsw) {
                auto& wsr = mswSpeciesRate->emplace_back(speciesName(sIdx));
                wsr.rate.reserve(eclWell.getConnections().size());
                for (std::size_t i = 0; i < eclWell.getConnections().size(); ++i) {
                    wsr.rate.emplace(eclWell.getConnections().get(i).segment(), 0.0);
                }
            }
        }

        // Get WSPECIES info
        std::vector<Scalar> wspecies(numSpecies());
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            wspecies[sIdx] = currentWSPECIES_(eclWell, speciesName(sIdx),
                                              simulator_.problem().wellModel().summaryState());
        }

        // Calculate well term
        const auto& ws = simulator_.problem().wellModel().wellState().well(well.name());
        for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
            // Get perforation rate
            const auto I = ws.perf_data.cell_index[i];
            const Scalar rate = well.volumetricSurfaceRateForConnection(I, waterPhaseIdx);

            // Injection
            if (rate > 0.0) {
                auto& source = getWellSource_(I);
                for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                    // Inject WSPECIES concentration
                    const Scalar inj_species_rate = rate * wspecies[sIdx];
                    source.speciesInjectionRate[sIdx] += inj_species_rate;

                    // Store for reporting here because WSPECIES is constant over time step
                    speciesRate[sIdx].rate += inj_species_rate;
                    if (isMsw) {
                        (*mswSpeciesRate)[sIdx].rate[eclWell.getConnections().get(i).segment()] += inj_species_rate;
                    }
                }
            }
            // Production
            // OBS: storing well rates for reporting done in transportSubstepExplicit_()
            else if (rate < 0.0) {
                // The species are produced with the concentration of the cell
                getWellSource_(I).waterProductionRate += rate;
                if (eclWell.isProducer()) {
                    producerConnections_.push_back(
                        {eclWell.seqIndex(),
                         I,
                         rate,
                         isMsw ? eclWell.getConnections().get(i).segment() : -1});
                }

                for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                    // Ensure reporting of cross-flow
                    const Scalar inj_species_rate = rate * wspecies[sIdx];
                    speciesRate[sIdx].rate += inj_species_rate;
                }
            }
        }
    }

    /*!
    * \brief Run geochemistry equilibrium solver
    *
    * \param dt Time step
    * \param I Cell index
    * \param initial_equil Optional bool indicating if this is an initial equilibrium solve
    * \param substepFraction Optional fraction of the time step at the end of the transport substep
    */
    void speciesEquationChemistryExplicit_(Scalar dt,
                                           unsigned I,
                                           bool initial_equil = false,
                                           Scalar substepFraction = 1.0)
    {
        // Ensure that dt == 0 for initial equilibration
        const double dt_geochem = initial_equil ? 0.0 : dt;

        // Initialize total and adsorbed concentrations before equilibration
        const auto& intQuants = simulator_.model().intensiveQuantities(I, 0);
        const auto& fs = intQuants.fluidState();
        double swat = decay<double>(fs.saturation(FluidSystem::waterPhaseIdx));

        // Water saturation at the end of the transport substep, from the same linear change of the
        // water volume over the time step as in the transport
        if (!initial_equil) {
            const Scalar volumeSubstep = vol1_[I] + substepFraction * (volumeNew_[I] - vol1_[I]);
            swat *= volumeSubstep / volumeNew_[I];
        }

        // Minimum water saturation in the initial equilibrium, so that the surface loading does not
        // depend on small saturations. The minerals change only with the water in the cell.
        const double swatMineral = swat;
        if (initial_equil) {
            swat = std::max<double>(swat, minSatInitial);
        }
        std::vector<double> Ctot;
        std::vector<double> Cads;
        const auto nAqu = geoChemInterface_->numberOfAqueousBasisSpecies();
        const auto nBasis = geoChemInterface_->numberOfBasisSpecies();
        Ctot.resize(nBasis, 0.0);
        Cads.resize(nBasis, 0.0);
        for (std::size_t k = 0; k < nAqu; ++k) {
            // Scale with water saturation, similar to ion exchangers
            if (!initial_equil) {
                Cads[k] = Cads_[k][I] / swat;
            }
            Ctot[k] = initial_equil ?
                concentrationInitial_[k][I] : concentration_[k][I];
            Ctot[k] += Cads[k];
        }

        // Setup mineral concentrations before equilibration
        // NOTE: assert checks that number of internal minerals are same as in the OPM deck
        const std::size_t nMin = geoChemInterface_->numberOfMinerals();
        assert(nMin == numMinerals());

        double* log_Amin_ptr = nullptr;
        double* Cmin_ptr = nullptr;
        std::vector<double> Cmin;
        const auto poro = decay<double>(intQuants.porosity());
        if (nMin > 0) {
            // Initial equilibrium solve requires calculation of initial mineral concentration
            if (initial_equil) {
                std::unordered_map<std::string, double> wt_frac;
                wt_frac.reserve(nMin);
                for (std::size_t l = 0; l < nMin; ++l) {
                    wt_frac.emplace(mineralNames_[l], minWt_[l][I]);
                }
                geoChemInterface_->calculate_initial_mineral_concentration(Cmin, poro, wt_frac);
                for (std::size_t l = 0; l < nMin; ++l) {
                    Cmin_[l][I] = Cmin[l];
                }
            }
            else {
                Cmin.resize(nMin);
                for (std::size_t l = 0; l < nMin; ++l) {
                    Cmin[l] = Cmin_[l][I];
                }
            }
            Cmin_ptr = Cmin.data();
            log_Amin_ptr = geoChemInterface_->get_log_a_mineral().data();
        }

        // No need for equilibrium solve in dry cells
        if (swat < tolSat) {
            return;
        }

        // Fluid properties
        const auto temp = decay<double>(fs.temperature(0));
        const auto pres = decay<double>(fs.pressure(waterPhaseIdx));

        // Oil and gas saturations
        double soil = 0.0;
        double sgas = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            soil = decay<double>(fs.saturation(FluidSystem::oilPhaseIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            sgas = decay<double>(fs.saturation(FluidSystem::gasPhaseIdx));
        }
        std::array<double, 3> mass_phase = { swat, soil, sgas };

        // Surface area
        double SA = geoChemInterface_->GetSurfaceArea();

        // Diffusion layer
        // OBS: calculated in set_surface_concentrations()!
        double frac_DL = 0.0;

        // Set surface concentration
        // NOTE: assert checks that number of internal ion exchange species are same as in the OPM deck
        const auto nIon = geoChemInterface_->numberOfIonExchange();
        assert(nIon == numIonEx());

        // Set concentrations for ion exchangers, surface complexes, and diffusive layers in Ctot
        std::unordered_map<std::string, double> Cion;
        if (nIon > 0) {
            Cion.reserve(nIon);
            for (std::size_t m = 0; m < nIon; ++m) {
                Cion.emplace(ionExNames_[m], Cio_[m][I]);
            }
        }
        geoChemInterface_->set_surface_concentrations(mass_phase[0], Ctot, frac_DL, Cion);

        // misc variables
        double pH = pH_[I];
        double sigma = sigma_[I];
        double psi = psi_[I];

        // Run equilibrium solver
        geoChemInterface_->SolveChem_I(Ctot.data(),
                                       Cads.data(),
                                       Cmin_ptr,
                                       log_Amin_ptr,
                                       temp,
                                       pres,
                                       poro,
                                       dt_geochem,
                                       SA,
                                       frac_DL,
                                       mass_phase,
                                       pH,
                                       sigma,
                                       psi);

        // Update mineral concentrations. The change is per water volume and Cmin_ is per pore
        // volume, so scale with water saturation
        if (nMin > 0) {
            for (std::size_t l = 0; l < nMin; ++l) {
                Cmin_[l][I] += Cmin[l] * swatMineral;
            }
        }

        // Update misc variables
        pH_[I] = pH;
        sigma_[I] = sigma;
        psi_[I] = psi;

        // Update species concentrations
        for (std::size_t k = 0; k < nAqu; ++k) {
            Ctot[k] -= Cads[k];
            Cads_[k][I] = Cads[k] * swat;

            if (initial_equil) {
                concentrationInitial_[k][I] = Ctot[k];
            }
            else {
                concentration_[k][I] = Ctot[k];
            }
        }
    }

private:
    std::vector<SpeciesVector> concentrationInitial_;
    std::vector<SpeciesVector> concentration_;
    std::vector<SpeciesVector> Cads_;
    std::vector<SpeciesVector> Cio_;
    std::vector<SpeciesVector> Cmin_;
    std::vector<SpeciesVector> minWt_;
    std::vector<double> pH_;
    std::vector<double> sigma_;
    std::vector<double> psi_;
    std::vector<Scalar> vol1_;
    bool initialEquilibrationDone_{false};

    std::vector<unsigned> interiorDofs_;
    std::vector<std::vector<FaceFlux>> faceFluxes_;
    std::vector<Scalar> volumeNew_;
    std::unordered_map<unsigned, WellSource> wellSources_;
    std::vector<ProducerConnection> producerConnections_;
    std::vector<SpeciesVector> concentrationNext_;
    std::vector<std::string> speciesNames_;
    std::vector<std::string> mineralNames_;
    std::vector<std::size_t> deckToSolverMineral_;
    std::vector<std::string> ionExNames_;

    std::shared_ptr<OpmGeoChemInterface> geoChemInterface_;

    std::unordered_map<int, std::vector<WellTracerRate<Scalar>>> wellSpeciesRate_;
    std::unordered_map<int, std::vector<MSWellTracerRate<Scalar>>> mSwSpeciesRate_;

    Simulator& simulator_;
    const EclipseState& eclState_;
    const CartesianIndexMapper& cartMapper_;
    ElementChunks<GridView, Dune::Partitions::All> element_chunks_;
    GeochemistryModelParameters<Scalar> param_;
};  // class GeochemistryModel
} // namespace Opm

#endif
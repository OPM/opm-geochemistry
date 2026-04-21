/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#include <opm/simulators/geochemistry/GeoChemIF.h>

#include <opm/simulators/geochemistry/Core/SimulationSchedule.hpp>

CGeoChemIF::CGeoChemIF(std::unique_ptr<InitChem> initChem)
    : geochemicalDatabaseReader_({}, false),
    ICS_full_(std::move(initChem))
{
    set_input_reader_keywords();
    if (ICS_full_)
    {
        allocate_memory_for_solver_and_splay_tree_etc();
    }
}

void CGeoChemIF::printCurrentVchemToFile(const std::string& file_name)
{
    GCS_->writeCurrentBasVecToFile(remove_suffixes(file_name) + "init.out");
}

void CGeoChemIF::printSplayTreeInfo() const
{
    solutions_manager_.print_information();
}

std::size_t CGeoChemIF::numberOfSolutions() const
{
    return ICS_sol_.size();
}

std::size_t CGeoChemIF::numberOfBasisSpecies() const
{
    return no_basis_species_;
}

std::size_t CGeoChemIF::numberOfAqueousBasisSpecies() const
{
    return no_aq_species_;
}

std::size_t CGeoChemIF::numberOfIonExchange() const
{
    return no_io_;
}

std::size_t CGeoChemIF::numberOfMinerals() const
{
    return no_minerals_;
}
int CGeoChemIF::getSolverCalls() const
{
    return GCS_->solver_state_.noCalls_;
}
bool CGeoChemIF::basis_species_added_by_geochem()
{
    if (ICS_full_->species_added_after_input_file_.size() > 0)
    {
        std::cout << "To be consistent with input minerals, add the following basis species to simulation:" << std::endl;
        for (auto& spec : ICS_full_->species_added_after_input_file_)
            std::cout << spec << std::endl;
        return true;

    }
    return false;
}

std::vector<std::string> CGeoChemIF::GetSpeciesNames(int component_type, const std::string& postfix) const
{
    auto add_postfix = [postfix](std::vector<std::string>& original_vector)
    {
      for(auto& element : original_vector) element += postfix;
    };

    auto names_of_species = ICS_full_->get_basis_species_of_type(component_type);
    if(!postfix.empty()) add_postfix(names_of_species);

    return names_of_species;
}

std::vector<std::string> CGeoChemIF::GetAdsorbedSpeciesNames() const
{
    // First ion exchangers, then surface complexes...
    auto all_adsorbed_species = this->GetSpeciesNames(GeochemicalComponentType::ION_EXCHANGE);
    auto surface_complexes = this->GetSpeciesNames(GeochemicalComponentType::SURFACE_COMPLEX);
    move_elements_from_vector_to_vector(surface_complexes, all_adsorbed_species);
    return all_adsorbed_species;
}

std::vector<std::string> CGeoChemIF::GetMineralNames() const
{
    return ICS_full_->mineral_name_;
}

std::vector<std::string> CGeoChemIF::GetBasisSpeciesNames() const
{
    return ICS_full_->basis_species_name_;
}

std::vector<std::string> CGeoChemIF::GetIONames() const
{
    return ICS_full_->io_name_;
}

std::vector<std::string> CGeoChemIF::GetSurfNames() const
{
    return ICS_full_->surf_name_;
}

double CGeoChemIF::GetSurfaceArea() const
{
    return ICS_full_->SA_;
}

double CGeoChemIF::GetDiffusionLayer() const
{
    return ICS_full_->d_DL_;
}

/**  If ICS_in is nullptr, use properties of ICS_full_. */
void CGeoChemIF::SetSurfaceConc(double swat,
                                std::vector<double>& C_tot,
                                double& frac_DL,
                                InitChem* ICS_in) const
{
    InitChem* ICS = ICS_in ? ICS_in : ICS_full_.get();

    const int firstIonExchangerIndex = ICS_full_->size_aq_;
    const double Cw_inv = 1.0 / swat;

    // Exchange species are located after aqueous species...
    for (int i = 0; i < ICS->size_io_; ++i)
    {
        C_tot[firstIonExchangerIndex + i] = ICS->c_io_[i] * Cw_inv;
    }
    // ... then surface complexes
    const double specific_surface_area = ICS->SA_;
    const double size_of_diffuse_layer = ICS->d_DL_;

    const int firstSurfaceComplexIndex = ICS_full_->size_aq_ + ICS_full_->size_io_;
    for (int i = 0; i < ICS->size_surf_; ++i)
    {
        C_tot[firstSurfaceComplexIndex + i] = ICS->c_surf_[i] * (1.0e18 / PhysicalConstants::Avogadro) * specific_surface_area * Cw_inv;
    }

    if (size_of_diffuse_layer > 0.0)
    {
        frac_DL = 1.0e-6 * specific_surface_area * size_of_diffuse_layer * Cw_inv;
    }
}

InitChemData* CGeoChemIF::getDataForSolution(const std::string& solutionName)
{
    for (auto& icsData : ICS_sol_data_)
    {
        if (icsData.solutionName_ == solutionName) return &icsData;
    }
    return nullptr;
}

/** Relative to ICS_full_. */
std::vector<std::vector<double>> CGeoChemIF::GetInitialSolutions() const
{
    std::vector<std::vector<double>> C_inj;
    C_inj.resize(ICS_sol_.size());

    for(std::size_t solIdx=0; solIdx < ICS_sol_.size(); ++solIdx)
    {
        // Make sure that species added by user are available
        C_inj[solIdx].resize(ICS_full_->size_aq_, 0.0);

        for (int specieIdx = 0; specieIdx < ICS_sol_[solIdx]->size_aq_; ++specieIdx)
        {
            int specieIdxFull = ICS_sol_[solIdx]->kcmap_[specieIdx];
            C_inj[solIdx][specieIdxFull] = ICS_sol_[solIdx]->c_vchem_[specieIdx];
        }
    }
    return C_inj;
}

void CGeoChemIF::initialize_geochemical_state(std::size_t no_grid_blocks)
{
    auto& state = grid_state_for_netcdf_;

    state.no_grid_blocks_ = no_grid_blocks;

    state.pH_.resize(no_grid_blocks, 7.0);
    state.psi_.resize(no_grid_blocks, 0.0);
    state.sigma_.resize(no_grid_blocks, 0.0);

    state.species_concentrations_.clear();
    for(const auto& bas_specie: ICS_full_->SM_basis_->row_name_)
    {
        state.species_concentrations_[bas_specie] = std::vector<double>(no_grid_blocks, 0.0);
    }
    for(const auto& sec_specie: ICS_full_->SM_all_->row_name_)
    {
        state.species_concentrations_[sec_specie] = std::vector<double>(no_grid_blocks, 0.0);
    }
}

/** It is the responsibility of the caller to ensure that block_index iscorrect, i.e., that it is the grid block currently being solved for. */
void CGeoChemIF::update_geochemical_state(std::size_t block_index)
{
    auto& state = grid_state_for_netcdf_;
    if(block_index >= state.no_grid_blocks_) return;

    const BasVec& Vchem = *GCS_->Vchem_;

    // Surface potential is calculated from concentration of basis specie E
    state.pH_[block_index] = Vchem.calc_pH();
    state.sigma_[block_index] = Vchem.calc_surface_charge();

    for(auto& [specie_name, block_concentrations]: state.species_concentrations_)
    {
        const double c_specie = Vchem.get_species_concentration(specie_name);
        block_concentrations[block_index] = c_specie;

        if(specie_name == "E")
        {
            const double RT_div_F = (PhysicalConstants::IdealGasConstant * Temp_) / PhysicalConstants::Faraday;
            if(c_specie <= 0.0) state.psi_[block_index] = 0.0;
            else state.psi_[block_index] = 1.0e3 * RT_div_F * std::log(c_specie); // mV
        }
    }
}

/**
 * TODO: The porosity parameter is currently unused.
 *
 * Solve geochemical problem over a time step.
 *
 * @param [in,out] ctot Total concentrations of basis species.
 * @param [in,out] ctot_ads Adsorbed concentrations.
 * @param [in,out] ctot_min Mineral concentrations.
 *                          NB: At exit, holds concentration changes (new-old).
 * @param log_a_min [in,out] Log activity of minerals (?).
 * @param Temp [in] Temperature in Kelvin.
 * @param Pres [in] Pressure in Pascal.
 * @param porosity [in] Should be a value in the range (0, 1].
 * @param dt [in] Time step in seconds.
 * @param SA [in] Specific surface area in m^2/L PV.
 * @param frac_DL [in] Volume fraction of diffuse layer.
 * @param mass_phase [in] Mass of water, oil, and gas phases.
 * @param pH_in [in,out] pH.
 * @param sigma_in [in,out] Surface charge density.
 * @param psi_in [in,out] Surface potential.
 */
void CGeoChemIF::SolveChem_I(double* ctot,
                             double* ctot_ads,
                             double* ctot_min,
                             double* log_a_min,
                             double Temp,
                             double Pres,
                             [[maybe_unused]] double porosity,
                             double dt,
                             double SA,
                             double frac_DL,
                             const std::array<double, 3>& mass_phase,
                             double& pH_in,
                             double& sigma_in,
                             double& psi_in)
{
    const bool INTERPOLATE  = (solutions_manager_.splayTreeResolution() > 0);
    const int num_val       = val_log_.size();

    const bool surface_species_are_present = ICS_full_->includes_ion_exchange()
                                             || ICS_full_->includes_surface_complexes();

    const int sizeb = GCS_->sizeB_;
    const int sizem = GCS_->sizeM_;

    GCS_->updateCurrentTimestep(dt);
    GCS_->diffuse_layer_props_.frac_DL_ = frac_DL;

    fill_zero(c_old_);
    fill_zero(val_log_);

    for (int i = 0; i < sizeb; ++i)
    {
        c_old_[i]                           = ctot[i];
    }

    for (int i = 0; i < sizem; ++i)
    {
        c_min_old_[i]                       = ctot_min[i];
        dc_min_[i]                          = 0.0;
    }

    BasVec** Vpp = solutions_manager_.get_vchem(ctot, ctot_min, ICS_full_.get());

    TreeNode* tn_tmp = nullptr;
    BasVec* Vp;

    int new_key = create_unique_key_from_positive_elements(ctot_min, sizem);
    int key;
    do
    {
        Vp = Vpp[new_key];
        if (dt == 0.0)
        {
            Vp = Vp->convert_nlin_rate_equations_to_equilibrium(ctot_min);
        }
        key = new_key;

        GCS_->Vchem_ = Vp;

        Vp->Temp_[gFluidPhase::WATER] = Temp;
        Vp->Pres_[gFluidPhase::WATER] = Pres;
        Vp->SA_ = SA;
        Vp->log_a_[Vp->pos_pH_] = -pH_in;
        Vp->log_m_[Vp->pos_pH_] = -pH_in;

        if (INTERPOLATE)
        {
            SplayTree& st_tmp = *solutions_manager_.get(key);

            std::vector<int> st_key = solutions_manager_.createKey(ctot, dt, Vp->Temp_[gFluidPhase::WATER], log_a_min);

            // Note: After the function call, st_key is in a moved-from state.
            const bool keyDidNotExistBefore = st_tmp.SplayNewKey(st_key);

            // After splaying, the new key is at the root. However, if the key
            // did NOT exist before, the root of the splay tree should now
            // point to the same address in memory as "new node".
            // Hence, we need to solve a new geochemical system...
            if (keyDidNotExistBefore)
            {
                const int new_mineral_key = GCS_->SolveChem_ST(
                    st_tmp.mineralCombinationKey(),
                    dt,
                    mass_phase,
                    ctot,
                    ctot_ads,
                    ctot_min,
                    log_a_min,
                    val_log_.data()
                );

                st_tmp.insertValues(new_mineral_key, val_log_.data());
            }
            tn_tmp = st_tmp.root();
            new_key = tn_tmp->mineral_combination_key_;
        }
        else
        {   // No interpolation
            new_key = GCS_->SolveChem_ST(
                key,
                0.0, // or dt...
                mass_phase,
                ctot,
                ctot_ads,
                ctot_min,
                log_a_min,
                val_log_.data()
            );

            if (!GCS_->massBalanceHasConverged())
            {
                for (int i = 0; i < sizeb; ++i)     ctot[i]     = c_old_[i];
                for (int i = 0; i < sizem; ++i)     ctot_min[i] = 0.0;
                return;
            }
        }

    } while (new_key != key);

    if (surface_species_are_present)
    {
        for (int i = 0; i < sizeb; ++i)
        {
            const double ctot_ads_old = ctot_ads[i];

            double dc;
            double dq;
            if(INTERPOLATE)
            {
                dc = tn_tmp->val_[i];
                dq = tn_tmp->val_[sizeb + i];
            }
            else
            {
                dc = val_log_[i];
                dq = val_log_[sizeb + i];
            }
            // TODO: What if ctot_ads[i] < 0 ?      [H+ can be negative]
            ctot_ads[i]     = ctot_ads_old - dq;
            ctot[i]         = c_old_[i] - dc;
            c_flux_[i]      = dc - dq;
        }
    }
    else  // No adsorption
    {
        for (int i = 0; i < sizeb; ++i)
        {
            const double dc = INTERPOLATE ? tn_tmp->val_[i] : val_log_[i];
            ctot[i]     = c_old_[i] - dc;
            c_flux_[i]  = dc;
        }
    }

    if (INTERPOLATE)
    {
        Vp->update_ctot_and_ctot_mineral(c_flux_, dc_min_);

        for (int i = 0; i < sizem; ++i)     ctot_min[i] = dc_min_[i];
    }
    else
    {
        for (int i = 0; i < sizem; ++i)     ctot_min[i] -= c_min_old_[i];
    }

    double* value_array = INTERPOLATE ? tn_tmp->val_.data() : val_log_.data();

    pH_in = value_array[num_val - 1];
    if (ICS_full_->includes_surface_complexes())
    {
        sigma_in    = value_array[2 * sizeb];
        psi_in      = value_array[2 * sizeb + 1];
    }

    if (ICS_full_->PRINT_DEBUG_CHEM_ > 0) {
        Vp->write_solution("Vp_sol.out");
    }
}

/** The first time it is called, allocates memory in ICS_sol_ for all POSSIBLE
 * solutions (may be less than the actual number).
 */
void CGeoChemIF::AddPhases(const PairOfPhaseNamesAndTypes& names_and_types, const std::string& unique_identifier)
{
    const std::size_t max_no_solutions = geoChemPhases_.phaseSize(GeochemicalPhaseType::AQUEOUS_SOLUTION);
    if (ICS_sol_.empty())
    {
        ICS_sol_.reserve(max_no_solutions);
    }
    else if(ICS_sol_.size() == ICS_sol_.capacity())
    {
        std::cout << "CGeoChemIF::AddPhases(): Cannot add InitChem object, need excess capacity in vector to avoid move constructor...\n";
        return;
    }
    else if(ICS_sol_.size() == max_no_solutions)
    {
        std::cout << "CGeoChemIF::AddPhases(): Cannot add InitChem object, have already read the maximum number of solutions...\n";
        return;
    }

    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const int max_size = geoChemPhases_.selectPhases(names_and_types, phases_to_be_used);

    // Construct temporary InitChem (unique) pointer, then immediately move it into the container
    auto ptr_to_new_ICS = InitChem::create_from_input_data(ICS_full_.get(), db_changes_made_by_user_, phases_to_be_used, max_size, unique_identifier);
    ICS_sol_.emplace_back(std::move(ptr_to_new_ICS));

    ICS_sol_data_.emplace_back(unique_identifier, no_aq_species_, no_minerals_);
}

void CGeoChemIF::AddAqueousSolution(const std::string& solution_name)
{
    const std::vector<std::string> phases_names = { solution_name };
    const std::vector<GeochemicalPhaseType> phases_types = { GeochemicalPhaseType::AQUEOUS_SOLUTION };

    AddPhases(std::make_pair(phases_names, phases_types), solution_name);  // use solution name as identifier for InitChem
}

void CGeoChemIF::set_input_reader_keywords()
{
    const std::vector<std::string> end_keywords = { R"(/end)", "end", R"(/ end)" };

    // (key, value)-pair keywords that are used by most solvers
    std::map<std::string, std::string> simple_key_value_pairs = std::map<std::string, std::string>();
    simple_key_value_pairs["add_species"] = R"("ions.txt")";
    simple_key_value_pairs["debug"] = "0";
    simple_key_value_pairs["interpolate"] = "0";
    simple_key_value_pairs["NumericalJacobian"] = "0";

    simple_key_value_pairs["force_mbal"] = "0";

    ptrInputReader_ = std::make_unique<InputReader>(simple_key_value_pairs, /*case_sensitive=*/ false);
    ptrInputReader_->define_block_keyword("chemtol", end_keywords);
    ptrInputReader_->define_block_keyword("inject", end_keywords);  // Only used by 1D transport solver

    // 4/7-22: Allow defining "BASIS_SPECIES", "SECONDARY_SPECIES", etc. directly in the main input file.
    for (const auto& keyword : GeochemicalDatabaseKeyword::ALL_DATABASE_KEYWORDS)
    {
        // When reading from top-level input file, we use this InputReader object:
        ptrInputReader_->define_block_keyword(keyword, end_keywords);
        // When reading from separate file (e.g., "ions.txt"), we use this InputReader object instead:
        geochemicalDatabaseReader_.define_block_keyword(keyword, end_keywords);
    }
}

/** Note that the interpolation flag must be set before calling this function.*/
void CGeoChemIF::allocate_memory_for_solver_and_splay_tree_etc()
{
    no_basis_species_   = ICS_full_->size_basis_;
    no_aq_species_      = ICS_full_->size_aq_;
    no_minerals_        = ICS_full_->size_min_;
    no_io_              = ICS_full_->size_io_;

    GCS_ = std::make_unique<GCSolver>(
        &logger_,
        nullptr,
        ICS_full_->size_basis_,
        ICS_full_->size_min_,
        ICS_full_->SM_basis_->noRows_,
        ICS_full_->SM_all_->noRows_,
        ICS_full_->SM_mineral_->noRows_
    );

    const std::size_t no_values_to_store = solutions_manager_.init(
        ICS_full_->INTERPOLATE_,
        ICS_full_->get_all_basis_species(),
        ICS_full_->get_all_minerals(),
        ICS_full_->includes_ion_exchange(),
        ICS_full_->includes_surface_complexes()
    );

    val_log_.resize(no_values_to_store, 0.0);

    c_old_.resize(no_basis_species_, 0.0);
    dc_min_.resize(no_basis_species_, 0.0);
    c_flux_.resize(no_basis_species_, 0.0);
    c_min_old_.resize(no_minerals_, 0.0);

    Cmin_.resize(no_minerals_, 0.0);
    Cchem_.resize(no_basis_species_, 0.0);
    Ccads_.resize(no_basis_species_, 0.0);
    Caq_.resize(no_basis_species_, 0.0);
    Cdl_.resize(no_basis_species_, 0.0);

    Charge_.resize(no_basis_species_, 0.0);
    for (int i = 0; i < no_basis_species_; ++i)
    {
        Charge_[i] = ICS_full_->SM_basis_->charge_[ICS_full_->pos_[i]];
    }

    mass_phase_ = std::array<double, 3>{1.0, 1.0, 1.0};
}

void CGeoChemIF::modify_selected_solver_options(GCSolver& GCS_in) const
{
    GCS_in.options_.CALC_JACOBI_NUM_ = std::stoi(ptrInputReader_->get_simple_keyword_value("NumericalJacobian"));

    const auto& [mbal_crit, pH_crit] = parse_chemtol(ptrInputReader_->get_block_keyword_content("chemtol"));
    GCS_in.options_.MBAL_CONV_CRITERION_ = mbal_crit;
    GCS_in.options_.PH_CONV_CRITERION_ = pH_crit;
}

void CGeoChemIF::modify_selected_solver_options_json(GCSolver& GCS_in, nlohmann::json& json_parsed) const
{
    if (json_parsed.contains("chemtol")) {
        const auto& [mbal_crit, pH_crit] = parse_chemtol(json_parsed["chemtol"]);
        GCS_in.options_.MBAL_CONV_CRITERION_ = mbal_crit;
        GCS_in.options_.PH_CONV_CRITERION_ = pH_crit;
    }
}

void CGeoChemIF::set_geochemical_database_modifications_from_user_input()
{
    //
    // Scans the previously read input for additional species, new chemical
    // reactions,  modifications to existing reactions, etc.
    //
    // Currently, there are two ways to do this:
    //
    //   (1) Add keywords like "BASIS_SPECIES", "SECONDARY_SPECIES", etc.
    //       directly to the main input file.
    //   (2) Add the keywords to a separate file, the filename of which is
    //       given by the keyword "ADD_SPECIES".
    //
    //                                      !!! IMPORTANT !!!
    //
    //   Option (1) takes precedence over option (2). In other words, if at
    //   least one of the keywords have been defined in the main input file,
    //   the "ADD_SPECIES" keyword is completely ignored.
    //

    bool has_modified_db = false;
    for (const auto& keyword : GeochemicalDatabaseKeyword::ALL_DATABASE_KEYWORDS)
    {
        const auto& key_content = ptrInputReader_->get_block_keyword_content(keyword);
        if (!key_content.empty())
        {
            db_changes_made_by_user_.insert_or_assign(keyword, key_content);
            has_modified_db = true;
        }
    }

    if (has_modified_db)
    {
        logger_.warning("Geochemical database was modified by keywords in top-level input file, ignoring contents of ADD_SPECIES...");
        return;
    }

    // IMPORTANT: We only use "ADD_SPECIES" if no keywords have been detected in the main input file...
    const std::string database_fn = remove_quotation_marks_from_filename(ptrInputReader_->get_simple_keyword_value("ADD_SPECIES"));
    std::ifstream db_file(database_fn, std::ios::binary);
    if (db_file.is_open())
    {
        geochemicalDatabaseReader_.read(db_file);
        for(const auto& keyword : GeochemicalDatabaseKeyword::ALL_DATABASE_KEYWORDS)
        {
            const auto& key_content = geochemicalDatabaseReader_.get_block_keyword_content(keyword);
            if (key_content.empty()) continue;
            db_changes_made_by_user_.insert_or_assign(keyword, key_content);
        }
    }
}

void CGeoChemIF::set_geochemical_database_modifications_from_json(const std::string& json_input)
{
    // Scans the previously read input for additional species, new chemical
    // reactions,  modifications to existing reactions, etc.
    //
    // Currently, there are two ways to do this:
    //
    //   (1) Add keywords like "BASIS_SPECIES", "SECONDARY_SPECIES", etc.
    //       directly to the main input file.
    //   (2) Add the keywords to a separate file, the filename of which is
    //       given by the keyword "ADD_SPECIES".
    //
    //                                      !!! IMPORTANT !!!
    //
    //   Option (1) takes precedence over option (2). In other words, if at
    //   least one of the keywords have been defined in the main input file,
    //   the "ADD_SPECIES" keyword is completely ignored.

    // Read JSON
    nlohmann::json json_parsed;
    if (!json_input.empty()) {
        std::ifstream in(json_input);
        if (in) {
            // Read JSON file
            in >> json_parsed;
        }
        else {
            // Parse JSON string
            json_parsed = nlohmann::json::parse(json_input);
        }
    }
    else {
        throw std::runtime_error("Error: No JSON input provided!");
    }

    // Insert database keywords from JSON input
    bool has_modified_db = false;
    for (const auto& keyword : GeochemicalDatabaseKeyword::ALL_DATABASE_KEYWORDS) {
        if (json_parsed.contains(keyword)) {
            const auto& content = json_parsed.at(keyword).get<std::vector<std::string>>();
            if (!content.empty()) {
                db_changes_made_by_user_.insert_or_assign(keyword, content);
                has_modified_db = true;
            }
        }
    }

    // Jump out early if database keywords have been encountered
    if (has_modified_db) {
        return;
    }

    // Read ADD_SPECIES file and add database keywords
    if (json_parsed.contains("ADD_SPECIES")) {
        const auto add_species_file = json_parsed.at("ADD_SPECIES").get<std::string>();
        std::ifstream db_file(add_species_file, std::ios::binary);
        if (db_file.is_open()) {
            geochemicalDatabaseReader_.read(db_file);
            for(const auto& keyword : GeochemicalDatabaseKeyword::ALL_DATABASE_KEYWORDS) {
                const auto& key_content = geochemicalDatabaseReader_.get_block_keyword_content(keyword);
                if (!key_content.empty()) {
                    db_changes_made_by_user_.insert_or_assign(keyword, key_content);
                }
            }
        }
    }

}

void CGeoChemIF::create_ICS_full(const std::vector<std::string>& species_names_in_order)
{
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const int max_size = geoChemPhases_.selectPhases(geoChemPhases_.getPhaseNamesAndTypes(), phases_to_be_used);

    ICS_full_ = InitChem::create_from_input_data
        (
            nullptr,
            db_changes_made_by_user_,
            phases_to_be_used,
            max_size,
            ICS_full_name_,
            species_names_in_order
        );
}

void CGeoChemIF::set_data_for_solution(const std::string& solutionName, double porosity)
{
    for (std::size_t solutionIndex = 0; solutionIndex < ICS_sol_.size(); ++solutionIndex)
    {
        if (ICS_sol_[solutionIndex]->name() == solutionName)
        {
            set_data_for_solution(solutionIndex, porosity);
            return;
        }
    }
    std::cout << "No solution named " << solutionName << "\n";
}

void CGeoChemIF::set_data_for_solution(std::size_t solIdx, double porosity)
{
    if (solIdx >= ICS_sol_.size()) return;  // throw error?

    const double r_density = ICS_full_->compute_rock_density(porosity);  // kg/m^3

    const InitChem& ICS = *ICS_sol_[solIdx];

    for (int i = 0; i < ICS.size_min_; ++i)
    {
        const auto name_of_buffer = ICS.get_mineral_name(i);
        const double input_concentration = ICS.c_mineral_[i];
        const double mol_weight = ICS.SM_mineral_->mol_weight_[i];

        double buffer_concentration = 0.0;
        if(is_gas_buffer(name_of_buffer))
        {
            buffer_concentration = input_concentration / mol_weight;
        }
        else
        {
            const double wt_fraction = input_concentration;
            const double fac = 1.0e-3*(r_density/mol_weight)*(1.0-porosity)/porosity;
            buffer_concentration = fac*wt_fraction;
        }

        const auto mineralIndex = ICS.krmap_[i];
        Cmin_[mineralIndex]                                 = buffer_concentration;
        ICS_sol_data_[solIdx].C_minerals_[mineralIndex]     = buffer_concentration;
        ICS_sol_data_[solIdx].Log_a_minerals_[mineralIndex] = ICS.log_a_mineral_[i];
    }

    for (int i = 0; i < ICS_sol_[solIdx]->size_aq_; ++i)
    {
        const int speciesIdx = ICS.kcmap_[i];
        ICS_sol_data_[solIdx].C_aq_[speciesIdx] = ICS.c_vchem_[i];
    }
}

/** Used for equilibrium calc. involving a single solution. */
void CGeoChemIF::equilibrate_solution(InitChem* ICS,
                                      bool equilibrate,
                                      const std::string& case_name,
                                      int GEOCHEMDEBUG,
                                      std::vector<double>& Ct,
                                      std::vector<double>& Cads,
                                      double& pH,
                                      std::vector<double>& scharge,
                                      BasVecInfo* basVecInfo,
                                      double frac_DL)
{
    if(!ICS)
    {
        std::cout << "CGeoChemIF::equilibrate_solution(): Need to input a non-null InitChem pointer!\n";
        return;
    }

    ICS->PRINT_DEBUG_CHEM_ = GEOCHEMDEBUG;

    if (!ICS->SM_basis_) throw InvalidInputException("CGeoChemIF::equilibrate_solution(): Database have not been correctly set, cannot find any basis species...");
    if (!case_name.empty() && GEOCHEMDEBUG) ICS->write(case_name + "_ICS.out");

    BasVec* Vp = BasVec::createFromInitChem(ICS);

    if (!case_name.empty() && GEOCHEMDEBUG) Vp->write(case_name + "_Vp.out");
   // removed when doing gas phase calculations
   // for (int i = 0; i < ICS->size_aq_; ++i)
   // {
   //     Vp->ctot_[ICS->pos_[i]] = Ct[ICS->kcmap_[i]];
   // }
   // Vp->add_gas_phase_to_total_conc();

    Vp->equilibrate_ = equilibrate;

    auto GCS = std::make_unique<GCSolver>(
        &logger_,
        Vp,
        Ct.size(),
        ICS->size_sup_min_ + ICS->size_rock_,
        ICS->SM_basis_->noRows_,
        ICS->SM_all_->noRows_,
        ICS->SM_mineral_->noRows_
    );

    // Note: In this function, we do not use either ICS_full_ or GCS_.
    //       Hence, we need to update the instances created in the current scope:
    modify_selected_solver_options(*GCS);
    ICS->possibly_change_inconsistent_options();

    if(frac_DL < 0)
    {   // Assume Sw=1
        GCS->diffuse_layer_props_.frac_DL_ = 1.0e-6 * Vp->SA_ * Vp->d_DL_;
    }
    else GCS->diffuse_layer_props_.frac_DL_ = frac_DL;

    bool converged = GCS->SolveChem();

    pH = -Vp->log_a_[Vp->pos_pH_];

    if (!scharge.empty())
    {
        scharge[0] = Vp->sch_;
        scharge[1] = Vp->psi_;
    }

    for (int i = 0; i < GCS->sizeB_; ++i)
    {
        Ct[i] = GCS->cT_[i];
        if (!Cads.empty())
        {
            Cads[i] = GCS->cT_ads_[i];
        }
    }

    if (!case_name.empty())
    {
        // Note: Here we use local GCS, not GCS_
        Vp->write_solution(case_name + "_solution.out", &GCS->diffuse_layer_props_);
        Vp->write_solution_chemistry(case_name + "_aq.out",/*include_reactions=*/true);
        Vp->write_buffers(case_name + "_buffer.out",/*include_reactions=*/true);
        Vp->write_tot_conc(case_name + "_species.out");
        if (Vp->size_gase_phase_ > 0)
        {
            Vp->write_gas_phase_solution(case_name + "_gas_phase.out");
        }

        GCS_->write_convergence_status("converged.txt", converged);
    }

    Vp->write_info(basVecInfo);

    if (GEOCHEMDEBUG >= 2 && !case_name.empty())
    {
        Vp->ICS_->SM_mineral_->write(case_name + "_mineral_table.out");
        Vp->ICS_->SM_all_->write(case_name + "_aq_table.out");
        Vp->ICS_->SM_basis_->write(case_name + "_basis_table.out");
    }

    //GCS->write_logK_to_file_T(298.15, 298.15 + 150, 5, 8e5, "XXtmp.out", GeochemicalComponentType::SURFACE_COMPLEX);

    delete Vp;
}

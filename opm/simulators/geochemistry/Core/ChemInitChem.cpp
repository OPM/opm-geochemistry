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
#include <opm/simulators/geochemistry/Core/ChemInitChem.h>

// CONSTRUCTORS

InitChem::InitChem(_constructor_tag,
    InitChem* ICS_full,
    const std::map<std::string, std::vector<std::string>>& db_changes_made_by_user,
    const std::map<GeochemicalPhaseType, GeoChemPhaseData>& phases_to_be_used,
    int max_size,
    const std::string& unique_identifier,
    const std::vector<std::string>& species_names_in_order,
    int interpolate)
{
    name_ = unique_identifier;
    INTERPOLATE_ = interpolate;

    // We have to do this in a very specific order...
    initialize_vectors(max_size);  // First, must allocate memory
    set_phase_composition(phases_to_be_used);
    reorder_species(species_names_in_order);

    if( db_changes_made_by_user.size() > 0)
    {
        std::cout << "Reading user changes to the default geochemical database...\n";
        read_geochemical_database(db_changes_made_by_user);
    }

    std::cout << "Initializing database..." << std::endl;
    init_chemistry_hkf();
    std::cout << "Initialization done!" << std::endl;

    // Need to update relative positions at the end (if applicable)
    if (ICS_full) set_relative_position_and_mineral_conc(*ICS_full);
}

InitChem::InitChem(
    int debug_mode,
    bool charge_balance,
    const ChemTable& SM_basis,
    const ChemTable& SM_all,
    const ChemTable& SM_mineral,
    const std::vector<std::string>& basis_name,
    const std::vector<double>& basis_conc,
    const std::vector<int>& basis_offset_index,
    const std::vector<std::string>& mineral_name,
    const std::vector<double>& mineral_conc,
    const std::vector<int>& mineral_offset_index,
    const std::vector<double>& mineral_log_af,
    const std::vector<int>& pos_gas,
    const std::vector<int>& phase,
    const std::vector<MineralRateLaw>& rate_laws,
    const std::vector<std::string>& io_name,
    const std::vector<std::string>& surf_name,
    const SurfaceChemistryOptions& surf_options,
    double specific_surface_area,
    double size_of_the_diffuse_layer
)
{
    PRINT_DEBUG_CHEM_ = debug_mode;
    charge_balance_ = charge_balance;

    // Probably larger than the actual size (we count some species twice)
    const int max_size = basis_name.size() +
        mineral_name.size() +
        io_name.size() +
        surf_name.size();

    initialize_vectors(max_size);

    size_basis_ = basis_name.size();
    size_io_ = io_name.size();
    size_surf_ = surf_name.size();
    size_aq_ = size_basis_ - size_io_ - size_surf_;
    assert(size_aq_ >= 1);

    size_min_ = mineral_name.size();
    size_sup_min_ = rate_laws.size();
    size_rock_ = size_min_ - size_sup_min_;
    size_gas_ = pos_gas.size();

    basis_species_name_ = basis_name;
    io_name_ = io_name;
    surf_name_ = surf_name;
    mineral_name_ = mineral_name;
    assert(mineral_name_.size() == static_cast<std::size_t>(size_min_));

    sup_min_name_.resize(size_sup_min_);
    for (int j = 0; j < size_sup_min_; ++j)
    {
        sup_min_name_[j] = mineral_name[j];
        rate_[j] = rate_laws[j];
    }

    buffer_name_.resize(size_rock_);
    for (int j = 0; j < size_rock_; ++j)
    {
        buffer_name_[j] = mineral_name[size_sup_min_ + j];
    }

    for (int i = 0; i < size_basis_; ++i)   c_vchem_[i] = basis_conc[i];
    for (int i = 0; i < size_basis_; ++i)   kcmap_[i] = basis_offset_index[i];
    for (int i = 0; i < size_min_; ++i)     krmap_[i] = mineral_offset_index[i];
    for (int i = 0; i < size_sup_min_; ++i) c_sup_min_[i] = mineral_conc[i];
    for (int i = 0; i < size_rock_; ++i)    c_buffer_[i] = mineral_conc[size_sup_min_ + i];
    for (int i = 0; i < size_rock_; ++i)    log_af_[i] = mineral_log_af[i];
    for (int i = 0; i < size_rock_; ++i)    phase_[i] = phase[i];
    for (int i = 0; i < size_gas_; ++i)     pos_gas_[i] = pos_gas[i];

    surface_options_from_input_ = surf_options;
    SA_ = specific_surface_area;
    d_DL_ = size_of_the_diffuse_layer;

    //  Note that we no longer use ICS_full as input here (directly)
    init_chemistry_hkf(SM_basis, SM_all, SM_mineral);
}

// end constructors

std::vector<std::string> InitChem::remove_organic_species()
{
    std::vector<std::string> species_to_remove;
    // Here we remove all species that contains e- != 0 and HCO3 != 0;
    int pos_e = SM_basis_->get_row_index("e-");
    int pos_hco = SM_basis_->get_row_index("HCO3-");
    if (pos_e > -1 && pos_hco > -1)
    {
        for (std::size_t i = 0; i < SM_all_->row_name_.size(); ++i)
        {
            if (SM_all_->M_[i][pos_e] != 0 && SM_all_->M_[i][pos_hco] != 0)
                species_to_remove.push_back(SM_all_->row_name_[i]);
        }
    }
    return species_to_remove;
}
InitChem* InitChem::create_from_ICS_full(const InitChem& ICS_full,
                                         const double* c_bas,
                                         const double* c_min,
                                         const std::vector<double>& cb_cp)
{
    // Add basis species present in the rock.
    std::vector<std::string> basis_name;
    std::vector<double> basis_conc;
    std::vector<int> basis_offset_index;

    basis_name.reserve(ICS_full.size_basis_);
    basis_conc.reserve(ICS_full.size_basis_);
    basis_offset_index.reserve(ICS_full.size_basis_);

    for (int i = 0; i < ICS_full.size_basis_; ++i)
    {
        if (cb_cp[i] != 0.0)
        {
            basis_name.push_back(ICS_full.basis_species_name_[i]);
            basis_conc.push_back(c_bas[i]);
            basis_offset_index.push_back(i);
        }
    }

    auto include_mineral = std::vector<bool>(ICS_full.size_min_, true);
    for (int i_min = 0; i_min < ICS_full.size_min_; ++i_min)
    {
        for (int i = 0; i < ICS_full.size_basis_; ++i)
        {
            const double stoichiometric_coefficient = ICS_full.SM_mineral_->M_[ICS_full.pos_min_[i_min]][ICS_full.pos_[i]];

            if (stoichiometric_coefficient != 0 && cb_cp[i] == 0)
            {   // Basis specie making up mineral does not exist!
                include_mineral[i_min] = false;
                break;
            }
        }
    }

    std::vector<std::string> mineral_name;
    std::vector<double> mineral_conc;
    std::vector<int> mineral_offset_index;
    std::vector<MineralRateLaw> rate_laws;
    std::vector<int> phase;
    std::vector<int> pos_gas;

    mineral_name.reserve(ICS_full.size_min_);
    mineral_conc.reserve(ICS_full.size_min_);
    mineral_offset_index.reserve(ICS_full.size_min_);
    rate_laws.reserve(ICS_full.size_sup_min_);
    phase.reserve(ICS_full.size_basis_);
    pos_gas.reserve(ICS_full.size_basis_);

    for (int i = 0; i < ICS_full.size_min_; ++i)
    {
        if(!include_mineral[i]) continue;

        // ...otherwise, we actually add the mineral
        mineral_name.push_back(ICS_full.mineral_name_[i]);
        mineral_conc.push_back(c_min[i]);

        if (i < ICS_full.size_sup_min_)
        {
            assert(ICS_full.rate_[i].has_value());
            rate_laws.emplace_back(ICS_full.rate_[i].value());
        }
        else
        {
            const int ph = ICS_full.phase_[i - ICS_full.size_sup_min_];
            phase.push_back(ph);

            if(ph == gFluidPhase::OIL || ph == gFluidPhase::GAS)
            {
                pos_gas.push_back(phase.size());
            }
        }
        mineral_offset_index.push_back(i);
    }

    // Specify which basis species are surface species
    std::vector<std::string> io_name;
    std::vector<std::string> surf_name;
    for (std::size_t i = 0; i < basis_name.size(); ++i)
    {
        const int type = ICS_full.SM_basis_->type_[ICS_full.pos_[basis_offset_index[i]]];

        if (type == GeochemicalComponentType::SURFACE_COMPLEX)      surf_name.push_back(basis_name[i]);
        else if (type == GeochemicalComponentType::ION_EXCHANGE)    io_name.push_back(basis_name[i]);
    }

    return new InitChem
        (
            ICS_full.PRINT_DEBUG_CHEM_,
            ICS_full.charge_balance_,
            *ICS_full.SM_basis_,
            *ICS_full.SM_all_,
            *ICS_full.SM_mineral_,
            basis_name,
            basis_conc,
            basis_offset_index,
            mineral_name,
            mineral_conc,
            mineral_offset_index,
            ICS_full.log_af_,
            pos_gas,
            phase,
            rate_laws,
            io_name,
            surf_name,
            ICS_full.surface_options_from_input_,
            ICS_full.SA_,
            ICS_full.d_DL_
        );
}


std::unique_ptr<InitChem> InitChem::create_from_input_data
    (
        InitChem* ICS_full,
        const std::map<std::string, std::vector<std::string>>& db_changes_made_by_user,
        const std::map<GeochemicalPhaseType, GeoChemPhaseData>& phases_to_be_used,
        int max_size,
        const std::string& unique_identifier,
        const std::vector<std::string>& species_names_in_order,
        int interpolate
    )
{
    return std::make_unique<InitChem>(_constructor_tag{},
                                      ICS_full,
                                      db_changes_made_by_user,
                                      phases_to_be_used,
                                      max_size,
                                      unique_identifier,
                                      species_names_in_order,
                                      interpolate);
}

void InitChem::possibly_change_inconsistent_options(bool enforce_charge_balance_when_using_diffuse_layer)
{
    if (surface_options_from_input_.ONLY_COUNTER_IONS)
    {
        surface_options_from_input_.INCLUDE_SURFACE_EXCESSES_IN_SURFACE_EQUATION = true;
    }

    if(enforce_charge_balance_when_using_diffuse_layer && surface_options_from_input_.CALCULATE_DIFFUSE_LAYER)
    {
        charge_balance_ = true;
        fmt::print("WARNING: Charge balance is enforced when we use the explicit diffusive layer model!\n.");
    }
}

int InitChem::read_geochemical_database(const std::map<std::string, std::vector<std::string>>& databaseReadFromInput)
{
    // Helper function to ensure that whenever a keyword is missing, we get back an empty vector.
    auto getKeywordContentPossiblyEmpty = [&databaseReadFromInput](const std::string& key)
    {
      if (databaseReadFromInput.find(key) == databaseReadFromInput.end()) return std::vector<std::string>();
      else return databaseReadFromInput.at(key);
    };

    basis_spec_db_ = getKeywordContentPossiblyEmpty(GeochemicalDatabaseKeyword::BASIS_SPECIES);
    sec_spec_db_ = getKeywordContentPossiblyEmpty(GeochemicalDatabaseKeyword::SECONDARY_SPECIES);
    surf_spec_db_ = getKeywordContentPossiblyEmpty(GeochemicalDatabaseKeyword::SURFACE_SPECIES);
    exch_spec_db_ = getKeywordContentPossiblyEmpty(GeochemicalDatabaseKeyword::EXCHANGE_SPECIES);
    min_phas_db_ = getKeywordContentPossiblyEmpty(GeochemicalDatabaseKeyword::MINERAL_PHASES);
    species_to_remove_ = getKeywordContentPossiblyEmpty(GeochemicalDatabaseKeyword::REMOVE_SPECIES);

    return 0;
}

void InitChem::set_relative_position_and_mineral_conc(InitChem& ICS)
{
    for (int c = 0; c < size_aq_; ++c)
    {
        // Aqueous species are stored first in the basis species vector
        kcmap_[c] = ICS.get_species_index(basis_species_name_[c]);
    }

    for (int c = 0; c < size_min_; ++c)
    {
        krmap_[c] = ICS.get_mineral_index(mineral_name_[c]);
        log_a_mineral_[c] = ICS.log_a_mineral_[c];
        c_mineral_[c] = ICS.c_mineral_[c];
    }
}

bool InitChem::includes_ion_exchange() const
{
    return (size_io_ > 0);
}

bool InitChem::includes_surface_complexes() const
{
    return (size_surf_ > 0);
}

double InitChem::compute_rock_density(double porosity) const
{
    bool error = false;

    // Along the way, checks input mineral weight fractions for consistency
    double tot_wt_frac_minerals = 0.0;
    double tot_volume_minerals = 0.0;

    for(int i_buf=0; i_buf < size_min_; ++i_buf)
    {
        const auto name_of_buffer = get_mineral_name(i_buf);

        // TODO: What about wt.-fraction for gases?
        if(is_gas_buffer(name_of_buffer))
        {
            continue;  // No contribution
        }

        const double wt_frac = c_mineral_[i_buf];
        if(wt_frac < 0 || wt_frac > 1.0)
        {
            fmt::print("Error: Weight fraction of mineral {} must be in the range [0, 1].\n", name_of_buffer);
            error = true;
        }
        if(wt_frac > 0 && porosity == 1.0)
        {
            fmt::print("Error: Cannot have minerals in the core when the porosity is equal to unity!\n");
            error = true;
        }

        if(wt_frac > 0.0)
        {
            const int buffer_idx = pos_min_[i_buf];
            assert(buffer_idx >= 0);
            const double mol_volume = SM_mineral_->mol_volume_[buffer_idx];
            const double mol_weight = SM_mineral_->mol_weight_[buffer_idx];

            tot_wt_frac_minerals += wt_frac;
            tot_volume_minerals += wt_frac*mol_volume/mol_weight;
        }
    }

    if(tot_wt_frac_minerals > 1.0)
    {
        fmt::print("Invalid weight-fraction of minerals in core, cannot sum to more than unity.\n");
        error = true;
    }
    if(error) std::exit(-1);

    static constexpr double average_rock_density = 2700.; // kg/m^3
    tot_volume_minerals += (1.0-tot_wt_frac_minerals) / average_rock_density;

    double rock_density = 0.0;
    if(tot_volume_minerals > NumericalConstants::ZERO_THRESHOLD)
    {
        rock_density = 1.0 / tot_volume_minerals;
    }
    return rock_density;
}

std::vector<std::string> InitChem::get_all_basis_species() const
{
    return basis_species_name_;
}

std::vector<std::string> InitChem::get_all_secondary_species() const
{
    return SM_all_->row_name_;
}

std::vector<std::string> InitChem::get_all_minerals() const
{
    return mineral_name_;
};

/**
 *  TODO: While most entries of basis_species_name_ are upper-case, surface complexes
 *        can be mixed-case...
 *
 * @return Not necessarily identical names as the corresponding rows in
 *         SM_basis_; the casing might be different...
 */
std::vector<std::string> InitChem::get_basis_species_of_type(int component_type) const
{
    std::vector<std::string> names;

    for(std::size_t i=0; i < basis_species_name_.size(); ++i)
    {
        // Note: SM_basis_ contains more species than basis_species_name_ (also H2O)
        if (SM_basis_->type_[pos_[i]] == component_type)
        {
            names.push_back(basis_species_name_[i]);
        }
    }
    return names;
}

std::vector<std::string> InitChem::get_secondary_species_of_type(int component_type) const
{
    const auto& all_names = SM_all_->row_name_;

    std::vector<std::string> names;
    for(std::size_t j=0; j < all_names.size(); ++j)
    {
        if(SM_all_->type_[j] == component_type)
        {
            names.push_back(SM_all_->row_name_[j]);
        }
    }
    return names;
}

std::string InitChem::name() const{ return name_; }

/**
 * Species can be entered as either Ca+2 or Ca.
 *
 * IMPORTANT NOTE (8/8-2022): The way the code is now, there is no guarantee that:
 *
 *      - Specie names will have consistent casing everywhere.
 *      - That specie names will have charges.
 *          --> E.g., we might end up comparing CA to Ca+2...
 *
 * At least for now, we avoid the issue by always comparing upper case
 * representations after first removing any charges.(not very efficient)
 *
 */
int InitChem::get_species_index(const std::string& specie_name) const
{
    const auto specie_name_without_charge = ParsedSpecieData(specie_name).specie_name_without_charge_;

    for (auto it = basis_species_name_.cbegin(); it != basis_species_name_.cend(); ++it)
    {
        const auto second_specie_without_charge = ParsedSpecieData(*it).specie_name_without_charge_;

        if (to_upper_case(second_specie_without_charge) == to_upper_case(specie_name_without_charge))
        {
            return static_cast<int>(it - basis_species_name_.cbegin());
        }
    }
    return -1;
}

int InitChem::get_mineral_index(const std::string& mineral_name) const
{
    // 8/8-22: To be on the safe side, always compare upper case names
    return index_of_string_in_vector_upper_case(mineral_name, mineral_name_);
}

std::string InitChem::get_basis_name(int pos) const
{
    if (!SM_basis_ || pos_.empty())
    {
        fmt::print("InitChem::get_basis_name(): No basis specie set at pos={:d}, returning...\n", pos);
        return {};
    }
    else
    {
        return SM_basis_->row_name_[pos_[pos]];
    }
}

std::string InitChem::get_mineral_name(int pos) const
{
    if (!SM_mineral_ || pos_min_.empty())
    {
        fmt::print("InitChem::get_mineral_name(): No mineral set at pos={:d}, returning...\n", pos);
        return {};
    }
    else
    {
        return SM_mineral_->row_name_[pos_min_[pos]];
    }
}

void InitChem::print_species_to_screen() const
{
    fmt::print("----Minerals----\n");
    for (int i = 0; i < size_min_; ++i)     fmt::print("{}\n", get_mineral_name(i));
    fmt::print("------------------------\n");

    fmt::print("----Basis species----\n");
    for (int i = 0; i < size_basis_; ++i)   fmt::print("{}\n", get_basis_name(i));
    fmt::print("------------------------\n");
}

void InitChem::write(const std::string& file_name) const
{
    std::ofstream fp(file_name, std::ios::out);

    const auto hyphens = std::string("------------------------------------------------------------------");

    fmt::print(fp, "Species Read:\n");
    for (int i = 0; i < size_basis_; ++i)   fmt::print(fp, "{:s}\t{:g}\n", basis_species_name_[i], c_vchem_[i]);
    fmt::print(fp, "\n{:s}\n\n", hyphens);

    fmt::print(fp, "No of mass balance species: {:d}\n", size_mass_);
    for (int i = 0; i < size_mass_; ++i)    fmt::print(fp, "{:s}\t", basis_species_name_[pos_mass_[i]]);
    fmt::print(fp, "\n{:s}\n\n", hyphens);

    fmt::print(fp, "No of equilibrium species: {:d}\nRockBuffer\n", size_rock_);
    for (int i = 0; i < size_rock_; ++i)    fmt::print(fp, "{:s}\n", buffer_name_[i]);
    fmt::print(fp, "\n{:s}\n\n", hyphens);

    fmt::print(fp, "No of non-linear minerals: {:d}\n", size_sup_min_);
    for (int i = 0; i < size_sup_min_; ++i)
    {
        const auto str_rate = convert_to_string(rate_[i].value());
        fmt::print(fp, "{:s} rate const:\t{:s}\n", sup_min_name_[i], str_rate);
    }
    fmt::print(fp, "\n{:s}\n\n", hyphens);
}

/** Aqueous species have to be read before surface species. */
void InitChem::initialize_vectors(int max_size)
{
    max_size_ = max_size;
    const int vector_size = max_size + 1;  // Why +1?

    pos_.resize(vector_size, -1);
    pos_mass_.resize(vector_size, -1);
    pos_buffer_.resize(vector_size, -1);
    pos_rock_.resize(vector_size, -1);
    pos_min_.resize(vector_size, -1);
    pos_sup_min_.resize(vector_size, -1);
    pos_min_bas_.resize(vector_size, -1);
    pos_sup_bas_.resize(vector_size, -1);
    pos_gas_.resize(vector_size, -1);
    pos_gas_phase_.resize(vector_size, -1);

    phase_.resize(vector_size, gFluidPhase::UNDEFINED);

    krmap_.resize(vector_size, -1);
    kcmap_.resize(vector_size, -1);

    c_vchem_.resize(vector_size, 0.0);
    c_ads_.resize(vector_size, 0.0);
    c_io_ =  std::vector<double>(vector_size, 0.0);
    c_surf_.resize(vector_size, 0.0);

    c_mineral_.resize(vector_size, 0.0);
    gas_phase_pressure_.resize(vector_size, 0.0);
    c_sup_min_.resize(vector_size, 0.0);
    c_buffer_.resize(vector_size, 0.0);
    log_a_mineral_.resize(vector_size, 0.0);
    log_af_sup_.resize(vector_size, 0.0);
    log_af_.resize(vector_size, 0.0);
    Sg_.resize(vector_size, 0.0);

    rate_.resize(vector_size);
}

void InitChem::reorder_species(const std::vector<std::string>& speciesNamesInOrder)
{
    if(vector_has_duplicates(to_lower_case(speciesNamesInOrder)))
    {
        fmt::print("InitChem::reorder_species(): Input vector of names has duplicates, do nothing...\n.");
        return;
    }

    std::size_t addedIndex = 0;
    for(const auto& specie: speciesNamesInOrder)
    {
        auto SPECIE = to_upper_case(specie); // basis_specie_name_ is always upper case
        const int swapIndex = index_of_element_in_container(SPECIE, basis_species_name_);
        if(swapIndex == -1) continue;
        std::swap(basis_species_name_[addedIndex], basis_species_name_[swapIndex]);
        std::swap(c_vchem_[addedIndex], c_vchem_[swapIndex]);
        ++addedIndex;
    }
}

void InitChem::set_phase_composition(const std::map<GeochemicalPhaseType, GeoChemPhaseData>& phases)
{
    // Finally, init from the previously stored data about solutions, equilibrium phases etc.
    for (const auto& [phaseType, phaseData] : phases)
    {
        for (const auto& [dummy, block_content] : phaseData)
        {
            switch (phaseType)
            {
            case GeochemicalPhaseType::AQUEOUS_SOLUTION:
            {
                read_solution(block_content);
                break;
            }
            case GeochemicalPhaseType::EQUILIBRIUM_MINERAL:
            {
                read_ieq(block_content);
                break;
            }
            case GeochemicalPhaseType::GAS_PHASE:
            {
                read_gas_phase(block_content);
                break;
            }
            case GeochemicalPhaseType::EXCHANGE_SITES:
            {
                read_io(block_content);
                break;
            }
            case GeochemicalPhaseType::SURFACE_COMPLEX:
            {
                read_surface(block_content);
                break;
            }
            case GeochemicalPhaseType::RATE_MINERAL:
            {
                read_rate(block_content);
                break;
            }
            }
        }
    }

    init_mineral_vector_in_correct_order();

}

/**
* Adds up supersaturated minerals (first) and equilibrium minerals / buffers (second).
* Assumes the user is lazy and have not input them in the required order.
*/
void InitChem::init_mineral_vector_in_correct_order()
{
    const auto no_minerals_in_input = sup_min_name_.size() + buffer_name_.size();

    if (mineral_name_.empty() && mineral_name_.size() != no_minerals_in_input)
    {
        int i_min = 0;
        for (const auto& mineral_name : sup_min_name_)
        {   // First, rate-minerals
            mineral_name_.push_back(mineral_name);
            log_a_mineral_[i_min] = log_af_sup_[i_min];
            c_mineral_[i_min] = c_sup_min_[i_min];
            ++i_min;
        }

        int i_buf = 0;
        for (const auto& buffer_name : buffer_name_)
        {   // Next, equilibrium
            mineral_name_.push_back(buffer_name);
            log_a_mineral_[i_min] = log_af_[i_buf];
            c_mineral_[i_min] = c_buffer_[i_buf];
            ++i_buf;
            ++i_min;
        }
        size_min_ = i_min;
        assert(size_min_ == static_cast<int>(mineral_name_.size()));
    }
}

void InitChem::init_chemistry_hkf()
{
    CP_ = ChemParam();

    Basis_db_   = std::make_unique<ChemTable>();
    Aq_db_      = std::make_unique<ChemTable>();
    Mineral_db_ = std::make_unique<ChemTable>();

    // Read databases (NB: Important to read basis species first)
    std::vector < std::pair < std::string, std::string >> bb = { std::make_pair("SiO2(AQ)","Si")};
    Basis_db_->initialize_full_database_hkf(ChemTable::Type::Basis,bb);
    Aq_db_->initialize_full_database_hkf(ChemTable::Type::Complex);
    Mineral_db_->initialize_full_database_hkf(ChemTable::Type::Mineral);
    // User-specified basis species
    std::vector<std::string> redox_decoupled_basis_species =
        {
            {"HSd- 4 33.0729 / HKF 2860 -3850	16.3	5.0119	4.9799	3.4765	-2.9849	3.42	-6.27	1.441 /"},
            {"Hdg 4 2.0158    /HKF 4236 -1000	13.8	5.1427	4.7758	3.8729	-2.9764	27.6251	5.093	-0.209/"},
            {"Ndg 4 28.0134  /HKF 4347 -2495	22.9	6.2046	7.3685	2.8539	-3.0836	35.7911	8.3726	-0.3468/"}
        };
    if (!basis_spec_db_.empty())
    {
        for (const auto& reaction: basis_spec_db_)    redox_decoupled_basis_species.push_back(reaction);

    }
    add_basis_species(redox_decoupled_basis_species, GeochemicalComponentType::AQUEOUS_COMPLEX);

    if (!surf_spec_db_.empty())
    {
        add_basis_species(surf_spec_db_,
                          GeochemicalComponentType::SURFACE_COMPLEX);
    }

    if (surface_options_from_input_.EXCHANGE_MODEL_DL)
    {
        exch_spec_db_.emplace_back("XDL-");
        exch_spec_db_.emplace_back("YDL+");
    }
    if (!exch_spec_db_.empty())
    {
        add_basis_species(exch_spec_db_,
                          GeochemicalComponentType::ION_EXCHANGE);
    }

    Basis_db_->init_tables();
    if (PRINT_DEBUG_CHEM_)
    {
        Basis_db_->write_csv("Basis_db.csv");
    }
    // User-specified secondary species
    std::vector<std::string> analytical_io_exch_spec =
        {
            {"NaX = X- + Na+ / ANA 0.0 /"},
            {"CaX2 = 2X- + Ca+2 / ANA -0.8 /"},
            {"MgX2 = 2X- + Mg+2 / ANA -0.6 /"},
            {"KX = X- + K+ / ANA -0.7 /"},
            {"BaX2 = 2X- + Ba+2 / ANA -0.91 /"},
            {"SrX2 = 2X- + Sr+2 / ANA -0.91 /"}
        };

    if (surface_options_from_input_.EXCHANGE_MODEL_DL)
    {
        // Cation exchange in the diffuse layer
        analytical_io_exch_spec.emplace_back("NaXDL = XDL- + Na+ / ANA 0.0 /");
        analytical_io_exch_spec.emplace_back("HXDL = XDL- + H+ / ANA 0.0 /");
        analytical_io_exch_spec.emplace_back("CaXDL2 = 2XDL- + Ca+2 / ANA -0.8 /");
        analytical_io_exch_spec.emplace_back("MgXDL2 = 2XDL- + Mg+2 / ANA -0.6 /");
        analytical_io_exch_spec.emplace_back("KXDL = XDL- + K+ / ANA -0.7 /");
        analytical_io_exch_spec.emplace_back("BaXDL2 = 2XDL- + Ba+2 / ANA -0.91 /");
        analytical_io_exch_spec.emplace_back("SrXDL2 = 2XDL- + Sr+2 / ANA -0.91 / ");

        // Anion exchange in the diffuse layer
        analytical_io_exch_spec.emplace_back("HCO3YDL = YDL+ + HCO3- / ANA 0.0 /");
        //analytical_io_exch_spec.push_back("CO3YDL2 = YDL+ + HCO3- - H+ / ANA 0.0 /");
        analytical_io_exch_spec.emplace_back("ClYDL = YDL+ + Cl- / ANA 0.0 /");
        analytical_io_exch_spec.emplace_back("SO4YDL2 = 2YDL+ + SO4-2 / ANA 0.0 /");
    }

    if (!sec_spec_db_.empty())
    {
        for (const auto& reaction: sec_spec_db_)    analytical_io_exch_spec.push_back(reaction);

        partition_vector_by_key(analytical_io_exch_spec, "HKF");  // Need to put HKF species first (!!!)
    }

    add_secondary_species(analytical_io_exch_spec);
    Aq_db_->init_tables();

    // User-specified minerals and gasses
    std::vector<std::string> redox_decoupled_species =
        {
            {"Hdg,g = Hdg	    /HKF 0	0	 0      31.234	6.52	0.78	0.12	3000/"},
            {"H2Sd,g = HSd- + H+ /HKF 0	-8021	-4931	49.185	7.81	2.96	-0.46	2300/"},
            {"Ndg,g = Ndg        /HKF 0	0	0	45.796	6.83	0.9	-0.12	3000/"}
        };
    if (!min_phas_db_.empty())
    {
        for (const auto& reaction: min_phas_db_)   redox_decoupled_species.push_back(reaction);
        partition_vector_by_key(min_phas_db_, "HKF");  // Need to put HKF species first (!!!)
    }
    add_mineral_phases(redox_decoupled_species);

    Mineral_db_->init_tables();
    Mineral_db_->add_cubic_EOS_parameters_for_gases();

    if (PRINT_DEBUG_CHEM_)
    {
        Aq_db_->write("Aq.out");
        Mineral_db_->write("Mineral.out");
        Basis_db_->write("Basis.out");
    }

    // ---> Add basis species
    std::vector<std::string> sim_bas_spec;

    for(const auto& specie: bspecies_always_included_)
    {
        sim_bas_spec.emplace_back(specie);
    }
    for(const auto& specie: basis_species_name_)
    {
        sim_bas_spec.emplace_back(specie);
    }



    // If a gas phase is defined
    check_gas_phase_for_equilibrium_phases(gas_phase_name_, buffer_name_);
    check_gas_phase_for_equilibrium_phases(gas_phase_name_, sup_min_name_);

    add_mineral_basis_species(buffer_name_, sim_bas_spec);
    add_mineral_basis_species(sup_min_name_, sim_bas_spec);
    add_mineral_basis_species(gas_phase_name_, sim_bas_spec);

    if (includes_ion_exchange())
    {
        add_surface_species(io_name_, sim_bas_spec);
    }

    if (includes_surface_complexes())
    {
        sim_bas_spec.emplace_back("E");
        add_surface_species(surf_name_, sim_bas_spec);
    }

    // To prevent species (e.g., H+) to be added multiple times.
    sim_bas_spec = remove_duplicates(sim_bas_spec);

    // Make reduced table (1 for basis species that are in the calculation, 0 otherwise)
    std::vector<int> basis_species_to_keep = Basis_db_->find_row_elements(sim_bas_spec);
    std::vector<int> columns_to_keep_for_basis_table(Basis_db_->noColumns_, 1);

    // All input given, now we do the transformation of basis species
    // basis_species_to_keep are unaltered as we put new basis species in the same place
    if (!basis_transformation_.empty())
        clean_transformation();
    if (!basis_transformation_.empty())
        basis_transformation();


    SM_all_ = Aq_db_->create_reduced_table(basis_species_to_keep, 1, ChemTable::skip_col_complex() + 1);
    SM_mineral_ = Mineral_db_->create_reduced_table(basis_species_to_keep, 1, ChemTable::skip_col_min() + 1);
    SM_basis_ = Basis_db_->create_reduced_table_rc(basis_species_to_keep, columns_to_keep_for_basis_table, 1, ChemTable::skip_col_bas() + 1);
    // remove some of the organic species ?


    SM_mineral_->remove_row_species(species_to_remove_);
    SM_basis_->remove_row_species(species_to_remove_);
    SM_all_->remove_row_species(species_to_remove_);

    if (REMOVE_ORGANIC_SPECIES)
    {
        std::vector<std::string> org_spec;
        org_spec = remove_organic_species();
        SM_all_->remove_row_species(org_spec);
    }

    init_mineral_vector_in_correct_order();

    if (!SM_basis_->mol_weight_.empty())
    {
        calculate_mol_weight_mineral(*SM_mineral_, *SM_basis_);
        // Also for the aqueous complexes (used for the solution density)
        calculate_mol_weight_mineral(*SM_all_, *SM_basis_);
    }

    if (SM_mineral_->type_.empty())
    {
        SM_mineral_->type_.resize(SM_mineral_->noRows_,
                                  GeochemicalComponentType::MINERAL);
    }

    SM_basis_->update_number_of_analytical_species();
    SM_all_->update_number_of_analytical_species();
    SM_mineral_->update_number_of_analytical_species();
    fix_hkf_units();

    SM_all_->set_up_sparse_matrix();
    SM_mineral_->set_up_sparse_matrix();
    SM_basis_->set_up_sparse_matrix();

    // Small consistency checks of the database
    check_db(SM_all_.get(), SM_basis_.get());
    check_db(SM_mineral_.get(), SM_basis_.get());

    // Check whether (nick)names are repeated.
    if (PRINT_DEBUG_CHEM_)
    {
        Basis_db_->check_for_duplicate_species(/* is_case_sensitive=*/ false, /* check_nick_names= */ false);
        Aq_db_->check_for_duplicate_species(/* is_case_sensitive=*/ false, /* check_nick_names= */ false);
        Mineral_db_->check_for_duplicate_species(/* is_case_sensitive=*/ false, /* check_nick_names= */ false);
        Basis_db_->check_for_duplicate_species(/* is_case_sensitive=*/ false, /* check_nick_names= */ true);
        Aq_db_->check_for_duplicate_species(/* is_case_sensitive=*/ false, /* check_nick_names= */ true);
        Mineral_db_->check_for_duplicate_species(/* is_case_sensitive=*/ false, /* check_nick_names= */ true);
    }

    // We only use the reduced basis in the rest of the program
    Basis_db_.reset();
    Aq_db_.reset();
    Mineral_db_.reset();

    finalize_init_chemistry_hkf();
}

void InitChem::init_chemistry_hkf(const ChemTable& SM_basis, const ChemTable& SM_all, const ChemTable& SM_mineral)
{
    // NOTE: Assumes some state have been set already...

    CP_ = ChemParam();

    std::vector<std::string> sim_bas_spec;
    for(const auto& specie: bspecies_always_included_)      sim_bas_spec.push_back(specie);
    for(const auto& specie: basis_species_name_)            sim_bas_spec.push_back(specie);
    // Note: If ion exchangers are included, they should already be present in basis_species_name_.
    if (includes_surface_complexes())
    {
        sim_bas_spec.emplace_back("E");
    }

    sim_bas_spec = remove_duplicates(sim_bas_spec);  // To prevent species (e.g., H+) to be added multiple times.

    // Make reduced table (1 for basis species that are in the calculation, 0 otherwise)
    std::vector<int> basis_species_to_keep = SM_basis.find_row_elements(sim_bas_spec);
    std::vector<int> columns_to_keep_for_basis_table(SM_basis.noColumns_, 1);

    SM_basis_ = SM_basis.create_reduced_table_rc(basis_species_to_keep, columns_to_keep_for_basis_table);
    SM_all_ = SM_all.create_reduced_table(basis_species_to_keep);
    SM_mineral_ = SM_mineral.create_reduced_table(basis_species_to_keep);

    SM_basis_->update_number_of_analytical_species();
    SM_all_->update_number_of_analytical_species();
    SM_mineral_->update_number_of_analytical_species();

    finalize_init_chemistry_hkf();

    if(PRINT_DEBUG_CHEM_)
    {
        write("test.out");
        SM_all_->write("a1.out");
        SM_mineral_->write("a2.out");
        SM_basis_->write("a5.out");
    }
}

void InitChem::finalize_init_chemistry_hkf()
{
    if (PRINT_DEBUG_CHEM_)
    {
        SM_basis_->write("SM_basis.out");
        SM_all_->write("SM_all.out");
        SM_mineral_->write("SM_mineral.out");
    }

    SM_all_->transpose();
    init_phase_mol_weights();
    chem_set_pos_relative_to_full_db();
    db_set_activity();
}

void InitChem::init_phase_mol_weights()
{
    //  - Correlation for oil: Mw = 470*gamma_oil^5.57, gamma_oil=0.8 "light oil".
    //  - For gas phase, use value for dry air.
    mol_weight_phase_[gFluidPhase::WATER] = SM_basis_->mol_weight_[SM_basis_->get_row_index("H2O")];
    mol_weight_phase_[gFluidPhase::OIL] = 0.470*std::pow(0.8, 5.57);
    mol_weight_phase_[gFluidPhase::GAS] = 28.9647e-3;
}

void InitChem::chem_set_pos_relative_to_full_db()
{
    for(int i=0; i < size_basis_; ++i)
    {
        const int index_in_basis = SM_basis_->get_row_index(basis_species_name_[i]);
        if(index_in_basis == -1)
        {
            error_and_exit("No basis specie named {:s} was found...", basis_species_name_[i]);
        }
        else pos_[i] = index_in_basis;
    }

    // Check if all the mineral species are valid species
    for(int i=0; i < size_sup_min_; ++i)
    {
        const int index_mineral = SM_mineral_->get_row_index(sup_min_name_[i]);
        if(index_mineral == -1)
        {
            error_and_exit("No rate-limited mineral named {:s} was found...", sup_min_name_[i]);
        }
        else
        {
            pos_sup_min_[i] = index_mineral;
            pos_min_[i]     = index_mineral;
        }
    }

    for (int i = 0; i < size_rock_; ++i)
    {
        const int index_buffer = SM_mineral_->get_row_index(buffer_name_[i]);
        if(index_buffer == -1)
        {
            error_and_exit("No buffer named {:s} was found...", buffer_name_[i]);
        }
        else
        {   // Note (15/1-2023): The same as in IORCoreSim
            pos_buffer_[i]              = index_buffer;
            pos_min_[size_sup_min_ + i] = index_buffer;
        }
    }

    for (int i = 0; i < size_gas_phase_; ++i)
    {
        const int index_gas = SM_mineral_->get_row_index(gas_phase_name_[i]);
        if (index_gas == -1)
        {
            error_and_exit("No gas phase {:s} was found...", gas_phase_name_[i]);
        }
        else
        {
            pos_gas_phase_[i] = index_gas;
            // Note (23/11-2024): @ah below is not correct, it should be possible to buffer same species
            // in both rock and gas, if gas volume is large mineral will be totally dissolved
            //pos_min_[size_sup_min_ + size_rock_ + i] = index_gas;
        }
    }
}

void InitChem::db_set_activity()
{
    if (SM_mineral_->log_af_.empty())       SM_mineral_->log_af_.resize(SM_mineral_->noRows_, 0.0);
    if (SM_mineral_->fugacity_.empty())     SM_mineral_->fugacity_.resize(SM_mineral_->noRows_, 1.0);
    for(int i=0; i < size_rock_; ++i)       SM_mineral_->log_af_[pos_buffer_[i]] = log_af_[i];
    for(int i=0; i < size_sup_min_; ++i)    SM_mineral_->log_af_[pos_sup_min_[i]] = log_af_sup_[i];


}

/**
* Units in HKF database are cal & bar - we use SI units.
* Must call this function after reading all species.
*/
void InitChem::fix_hkf_units() const
{
    SM_basis_->fix_hkf_units(ChemTable::Type::Basis);
    SM_all_->fix_hkf_units(ChemTable::Type::Complex);
    SM_mineral_->fix_hkf_units(ChemTable::Type::Mineral);
}

void InitChem::add_basis_species(const std::vector<std::string>& buffer, int type)
{
    if (Basis_db_->noColumns_ == 0) // database is empty
    {
        Basis_db_->resize(0, 15);  // Bugfix (20/8-22): 14 --> 15.
        Basis_db_->col_name_[0].assign("type");
        Basis_db_->col_name_[1].assign("charge");
        Basis_db_->col_name_[2].assign("scharge");
        Basis_db_->col_name_[3].assign("mol_weight");
        Basis_db_->col_name_[4].assign("a0");
        Basis_db_->col_name_[5].assign("DeltaG");
        Basis_db_->col_name_[6].assign("DeltaH");
        Basis_db_->col_name_[7].assign("S");
        Basis_db_->col_name_[8].assign("a1");
        Basis_db_->col_name_[9].assign("a2");
        Basis_db_->col_name_[10].assign("a3");
        Basis_db_->col_name_[11].assign("a4");
        Basis_db_->col_name_[12].assign("c1");
        Basis_db_->col_name_[13].assign("c2");
        Basis_db_->col_name_[14].assign("omega");
    }

    std::vector<std::string> basis;
    basis.reserve(buffer.size());

    for (const auto& str_basis_specie : buffer)
    {
        double a0 = (type == GeochemicalComponentType::AQUEOUS_COMPLEX) ? DEFAULT_ION_SIZE : 0.0;
        double mol_weight = 0.0;
        double charge = 0.0;
        double surf_charge = 0.0;

        std::stringstream iss(substring_between(str_basis_specie, "", "/"));
        std::string bn;
        iss >> bn >> a0 >> mol_weight;

        const auto specieData = ParsedSpecieData(bn);
        if (type == GeochemicalComponentType::AQUEOUS_COMPLEX || type == GeochemicalComponentType::ION_EXCHANGE){
            charge = specieData.charge_;
        }
        else if (type == GeochemicalComponentType::SURFACE_COMPLEX){
            surf_charge = specieData.charge_;
        }

        std::vector<double> row_val;
        row_val.resize(5);  // at least these 5:
        row_val[0] = type;
        row_val[1] = charge;
        row_val[2] = surf_charge;
        row_val[3] = mol_weight;
        row_val[4] = a0;

        basis.emplace_back(specieData.specie_name_);

        const std::size_t indexHKF = str_basis_specie.find("HKF");
        const std::size_t indexANA = str_basis_specie.find("ANA");

        LogKModel model = LogKModel::HKF;
        if (indexHKF != std::string::npos)
        {
            model = LogKModel::HKF;
            const auto str_values = substring_between(str_basis_specie, "HKF", "/");
            for (const auto& value: parse_doubles(str_values)) row_val.push_back(value);
        }
        else if (indexANA != std::string::npos)
        {
            model = LogKModel::ANA;
            row_val.push_back(0.0);
        }
        Basis_db_->add_species_row(specieData.specie_name_, row_val);
        Basis_db_->model_[Basis_db_->noRows_ - 1] = model;
    }
    // Added for cases where user enters same basis specie multiple times
    std::set<int> uniqueBasis;
    for(const auto& name: basis)
    {
        uniqueBasis.insert(Basis_db_->get_row_index(name));
    }
    std::vector<std::string> basis_unique;
    for (int num : uniqueBasis)
    {
        basis_unique.push_back(Basis_db_->row_name_[num]);
    }
    Aq_db_->add_col(basis_unique);
    Mineral_db_->add_col(basis_unique);
}

/** NB: Needs to be called after add_basis_species().*/
void InitChem::add_secondary_species(const std::vector<std::string>& buffer)
{
    for (const auto& str_secondary_specie: buffer)
    {
        auto [vector_of_species_data, logK_data] = parseSecondarySpecie(str_secondary_specie);
        const auto& secondary_specie = vector_of_species_data[0];

        std::vector<double> charge = { secondary_specie.charge_ };
        std::vector<std::string> list_of_basis_species;
        std::vector<double> stoichiometric_coefficients;

        for (std::size_t basisIndex=1; basisIndex < vector_of_species_data.size(); ++basisIndex)
        {
            const auto& basis_specie = vector_of_species_data[basisIndex];
            list_of_basis_species.push_back(basis_specie.specie_name_);
            charge.push_back(basis_specie.charge_);
            stoichiometric_coefficients.push_back(basis_specie.coefficient_);
        }
        assert(charge.size() == list_of_basis_species.size() + 1);
        assert(stoichiometric_coefficients.size() ==  list_of_basis_species.size());

        // The rest of the code is used to check that basis species are ok + extract values for surface and exchange species
        bool all_basis_species_are_valid = true;
        bool add_surface_complex_species = false;
        bool add_ion_exchange_species = false;

        double ion_exchange_charge = 0.0;
        double a0_surface = 0.0;
        double E_coeff = 0.0;  // to count the coefficient of E0 (if a surface complexation reaction)

        for (std::size_t basisIndex=1; basisIndex < vector_of_species_data.size(); ++basisIndex)
        {
            const auto& basis_specie = vector_of_species_data[basisIndex];

            const int pos = Basis_db_->get_row_index(basis_specie.specie_name_);
            if (pos == -1)
            {
                all_basis_species_are_valid = false;
                if (PRINT_DEBUG_CHEM_ >= DebugInfoLevel::ALOT)
                {
                    fmt::print(
                        "geochem::add_species(): I am not adding secondary species {}, because basis species {} does not exist in the database!",
                        secondary_specie.specie_name_,
                        basis_specie.specie_name_
                    );
                }
            }
            else
            {
                const auto& type = Basis_db_->type_[pos];
                if (type == GeochemicalComponentType::AQUEOUS_COMPLEX)
                {
                    a0_surface = Basis_db_->a0_[pos];
                    ion_exchange_charge = Basis_db_->charge_[pos];  // NB: Assumes there is only one ion in the (half-)reaction
                    E_coeff -= basis_specie.charge_;
                }
                if (type == GeochemicalComponentType::SURFACE_COMPLEX)
                {
                    add_surface_complex_species = true;
                }
                if (type == GeochemicalComponentType::ION_EXCHANGE)
                {
                    add_ion_exchange_species = true;
                }
            }
        } // end for (basis loop)

        if (all_basis_species_are_valid)
        {
            int type;
            double a0;
            double charge_complex;

            if(add_surface_complex_species)
            {
                type = GeochemicalComponentType::SURFACE_COMPLEX;
                a0 = 0.0;
                charge_complex = vector_of_species_data[0].charge_;

                list_of_basis_species.emplace_back("E");
                stoichiometric_coefficients.push_back(E_coeff);
            }
            else if(add_ion_exchange_species)
            {
                type = GeochemicalComponentType::ION_EXCHANGE;
                a0 = a0_surface;  // E.g., if NaX = Na+ + X-, use a0 for Na+
                charge_complex = ion_exchange_charge;
            }
            else
            {
                type = GeochemicalComponentType::AQUEOUS_COMPLEX;
                a0 = DEFAULT_ION_SIZE;
                charge_complex = vector_of_species_data[0].charge_;
            }

            Aq_db_->add_species_row(secondary_specie.specie_name_, {});
            Aq_db_->add_reaction(
                secondary_specie.specie_name_,
                type,
                logK_data.model_,
                charge_complex,
                a0,
                /*mol_volume=*/ 0.0,
                logK_data.values_,
                list_of_basis_species,
                stoichiometric_coefficients);

        }  // end isValid()
    }  // end current secondary specie
}

/** Note: If a new species is detected, adds it to species_name. */
void InitChem::add_surface_species(const std::vector<std::string>& buffer,
                                   std::vector<std::string>& bnames)
{
    for (const auto& buffer_name: buffer)
    {
        const bool found = element_is_in_container(buffer_name,
                                                   basis_species_name_);

        if (!found)
        {
            basis_species_name_.push_back(buffer_name);
            bnames.push_back(to_upper_case(buffer_name));
            pos_mass_[size_mass_] = size_basis_;
            ++size_basis_;
            ++size_mass_;
        }
    }
}

/**
 * NB: Needs to be called after add_basis_species().
 *
 * HKF parameters:  [mol_volume, DeltaG, DeltaH, S, a1, a2, a3, a7].
 * ANA:             [mol_volume, a0, a1, a2, a3, a4], where logK =  a0 + a1*T + a2/T + a3*logT + a4/T^2.
*/
void InitChem::add_mineral_phases(const std::vector<std::string>& buffer)
{
    for (const auto& str_mineral_phase : buffer)
    {
        auto [vector_of_species_data, logK_data] = parseMineralPhase(str_mineral_phase);
        const auto& mineral_phase = vector_of_species_data[0];

        bool all_basis_species_are_valid = true;
        std::vector<std::string> list_of_basis_species;
        std::vector<double> stoichiometric_coefficients;

        for (std::size_t basisIndex = 1; basisIndex < vector_of_species_data.size(); ++basisIndex)
        {
            const auto& basis_specie = vector_of_species_data[basisIndex];
            const int pos = Basis_db_->get_row_index(basis_specie.specie_name_);
            if (pos == -1)
            {
                all_basis_species_are_valid = false;
                fmt::print("geochem::add_mineral_phases(): I am not adding mineral {:s}, because basis species {:s} does not exist in the database!\n",
                           mineral_phase.specie_name_,
                           basis_specie.specie_name_);
            }
            else
            {
                list_of_basis_species.push_back(basis_specie.specie_name_);
                stoichiometric_coefficients.push_back(basis_specie.coefficient_);
            }
        }

        if (all_basis_species_are_valid)
        {
            Mineral_db_->add_species_row(mineral_phase.specie_name_, {});
            Mineral_db_->add_reaction(
                mineral_phase.specie_name_,
                0, // NB: Should ideally be GeochemicalComponentType::MINERAL, but not sure that works...
                logK_data.model_,
                0.0, // charge
                DEFAULT_ION_SIZE,
                logK_data.mol_volume_,
                logK_data.values_,
                list_of_basis_species,
                stoichiometric_coefficients);
        }

    } // end current mineral phase
}
void InitChem::check_gas_phase_for_equilibrium_phases(const std::vector<std::string>& gas_phase, const std::vector<std::string>& buffer)
{
    if (gas_phase.size() < 1)
        return;
    for (const auto& gas_nm : gas_phase)
    {
        const int index_of_gas = Mineral_db_->get_row_index(gas_nm);
        if (index_of_gas < 0)
        {
            error_and_exit("Gas phase {:s}, does not exist in the geochemical database!", gas_nm);
        }
        for (const auto& buf_nm : buffer)
        {
            const int index_of_mineral = Mineral_db_->get_row_index(buf_nm);
            if (index_of_gas == index_of_mineral)
            {
                error_and_exit("Input error: {:s} cannot be defined as mineral phase and at same time defined in the gas phase!", gas_nm);
            }
        }
    }
}


int InitChem::add_mineral_basis_species(const std::vector<std::string>& mineral_names, std::vector<std::string>& bnames)
{
    int no_added_species = 0;

    if(!Mineral_db_)
    {   // Mineral database has not been initialized yet, do nothing...
        return no_added_species;
    }

    for (const auto & mineral_name : mineral_names)
    {
        const int index_of_mineral = Mineral_db_->get_row_index(mineral_name);
        if(index_of_mineral < 0)
        {
            error_and_exit("Tried to add mineral basis species for {:s}, but the mineral does not exist in the geochemical database!", mineral_name);
        }

        // Basis species starts later, see buffer database...
        for (int j = ChemTable::skip_col_min(); j < Mineral_db_->noColumns_; ++j)
        {
            if (fabs(Mineral_db_->M_[index_of_mineral][j]) > 0)
            {
                if (!basis_specie_defined(j - ChemTable::skip_col_min(), bnames))
                {
                    // NOTE: At the moment this may be printed out multiple times (for different InitChem objects)
                    if (PRINT_DEBUG_CHEM_ >= DebugInfoLevel::ALOT)
                    {
                        fmt::print(
                            "You have added mineral {}, but basis species {} is missing. Trying to add basis specie {} to solution...",
                            mineral_name,
                            Basis_db_->row_name_[j - ChemTable::skip_col_min()],
                            Basis_db_->row_name_[j - ChemTable::skip_col_min()]
                        );
                    }

                    // 8/8-22: Be more consistent about using upper case.
                    const auto& name_of_added_basis_specie = to_upper_case(Basis_db_->row_name_[j - ChemTable::skip_col_min()]);
                    bnames.push_back(name_of_added_basis_specie);
					species_added_after_input_file_.push_back(name_of_added_basis_specie);
                    basis_species_name_.push_back(name_of_added_basis_specie);
                    c_vchem_[size_aq_] = DefaultValues::MolarConcentration;
                    ++size_aq_;
                    ++size_basis_;
                    ++no_added_species;
                }
            }
        }
    }
    return no_added_species;
}

std::vector<std::string> InitChem::read_specie(const std::string& basis_specie_name,
                                  const std::string& line,
                                  double& conc) const
{
    const std::vector<std::string> lineTab = split_string(line);
    if (lineTab.empty())
    {
        conc = DefaultValues::MolarConcentration;
        if (PRINT_DEBUG_CHEM_)
        {
            fmt::print("No concentration given for {}, set to {} mol/L.\n",
                       basis_specie_name,
                       conc);
        }
    }
    else {
        conc = std::stod(lineTab[0]);
    }

    if(lineTab.size() >1)
    {
        std::vector<std::string> lt = std::vector<std::string>(lineTab.begin() + 1, lineTab.end());
        return to_upper_case(lt);
    }
    return std::vector<std::string> {""};
}

int InitChem::read_solution(const GeoChemBlockContent& block_content)
{
    if (block_content.empty())
        return 0;
    else
    {
        double current_basis_concentration = 0.0;
        int current_basis_index = size_basis_;
        for (const auto& [input_name, input_value] : block_content)
        {
            const std::string basis_specie_name = to_upper_case(input_name);
            const std::vector<std::string> w = read_specie(basis_specie_name,
                                              input_value,
                                              current_basis_concentration);

         //   const int pos = (basis_specie_name == "PH")
         //                   ? in_list("H", basis_species_name_)
         //                   : in_list(basis_specie_name, basis_species_name_);
            int pos;
            if (basis_specie_name == "PH")
            {
                pos = in_list("H", basis_species_name_);
            }
            else if (basis_specie_name == "PE")
            {
                pos = in_list("E-", basis_species_name_);
                current_basis_concentration = POW10(-current_basis_concentration);
            }
            else
                pos = in_list(basis_specie_name, basis_species_name_);

            // Ignore species with a concentration that is (practically) zero?
            const bool species_not_added = (pos == -1);

            if (species_not_added
                && current_basis_concentration > POW10(NumericalConstants::LOG10_ALMOST_ZERO))
            {
                std::string species_to_add= basis_specie_name;
                if (basis_specie_name == "PH")
                {
                    species_to_add = "H";
                }
                else if(basis_specie_name == "PE")
                {
                    species_to_add = "E-";
                };

                basis_species_name_.push_back(species_to_add);

                if (basis_specie_name == "PH")
                {
                    if (current_basis_concentration < 1)
                    {   // Unrealistic pH value... => Use default
                        current_basis_concentration = DefaultValues::pH;
                    }
                    current_basis_concentration = pow(10.0, -current_basis_concentration) + WHTO_;
                }
 //               else if (basis_specie_name == "PE")
 //               {
 //                   if (current_basis_concentration < 1) Maybe change in future??
 //                   {   // Unrealistic pe value... => Use default
 //                       current_basis_concentration = DefaultValues::pe;
 //                   }
 //                   current_basis_concentration = pow(10.0, -current_basis_concentration);
 //               }

                c_vchem_[current_basis_index] = current_basis_concentration;

                // Note: We have to increment the basis specie index at the end
                ++size_basis_;
                ++size_aq_;
                ++current_basis_index;

                // If at least one of the solutions has the keyword "CHARGE"
                if (w[0] == "CHARGE") charge_balance_ = true;

                if (w[0] == "FIX")
                {
                    keep_ph_fixed_ = true;
                }

                if (w[0] == "AS") // trigger basis transformation
                {
                    if (w.size() > 1)
                    {
                        basis_transformation_[species_to_add] = w[1];
                    }
                    else
                    {
                        fmt::print("Keyword \"as\" needs to acompanied by basis name\n");
                        fmt::print("Usage: old basis name --> new basis name\n");
                    }
                }
            }
        }
        return size_basis_;
    }
}

int InitChem::read_rate(const GeoChemBlockContent& block_content)
{
    if (block_content.empty())
        return 0;
    else
    {
        const int i0 = size_sup_min_; // some might have been read earlier

        for(const auto &[mineral_name, mineral_data]: block_content)
        {
            if (in_list(mineral_name, sup_min_name_) == -1)
            {
                c_sup_min_[size_sup_min_] = 0.0;  // Default wt-fraction is zero

                // Set default rate parameters (in order Sg, log_af, logEa_1, k_1, logEa_2, k_2, m, n, n_acid)
                // REFACTORING POTENTIAL: Here we might potentially override the defaults set in MineralRateLaw...
                std::vector<double> rate_params_current_mineral = std::vector<double>(9, 0.0);
                rate_params_current_mineral[0] = 1.0;  // Sg
                rate_params_current_mineral[3] = 1.0e5;  // k_1
                rate_params_current_mineral[6] = 1.0;  // m
                rate_params_current_mineral[7] = 1.0; // m
                rate_params_current_mineral[8] = 1.0; // n_acid

                const std::vector<std::string> lineTab = split_string(mineral_data);

                if (!lineTab.empty())
                {
                    const auto mineral_wt_fraction = std::stod(lineTab[0]);
                    c_sup_min_[size_sup_min_] = mineral_wt_fraction;
                }

                sup_min_name_.push_back(mineral_name);

                if(lineTab.size() > 1)
                {
                    int parameterIndex = 0;
                    for (auto it = lineTab.cbegin() + 1; it < lineTab.cend(); ++it)
                    {
                        rate_params_current_mineral[parameterIndex] = std::stod(*it);
                        ++parameterIndex;
                    }

                }

                rate_[size_sup_min_].emplace();  // constructs a MineralRateLaw object with default parameters

                Sg_[size_sup_min_] = rate_params_current_mineral[0]; // Sg
                rate_[size_sup_min_]->SA_ = rate_params_current_mineral[0]; // Sg
                log_af_sup_[size_sup_min_] = rate_params_current_mineral[1];  // log_af
                rate_[size_sup_min_]->k1_ = rate_params_current_mineral[3];
                rate_[size_sup_min_]->k2_ = rate_params_current_mineral[5];
                rate_[size_sup_min_]->m_ = rate_params_current_mineral[6];
                rate_[size_sup_min_]->n_ = rate_params_current_mineral[7];
                rate_[size_sup_min_]->Ea1_ = rate_params_current_mineral[2];
                rate_[size_sup_min_]->Ea2_ = rate_params_current_mineral[4];
                rate_[size_sup_min_]->n_acid_ = rate_params_current_mineral[8];
                ++size_sup_min_;
            }
        }
        return size_sup_min_ - i0;
    }
}

int InitChem::read_ieq(const GeoChemBlockContent& block_content)
{
    if (block_content.empty())
        return 0;
    else
    {
        int i0 = size_rock_; // some might have been read in solution keyword
        for (const auto& [buf_name, buf_values] : block_content)
        {
            if (in_list(buf_name, buffer_name_) == -1) // not added before
            {
                const std::vector<std::string> lineTab = split_string(buf_values);

                buffer_name_.push_back(buf_name);
                phase_[size_rock_] = gFluidPhase::WATER;
                if (lineTab.empty())
                {
                    if (PRINT_DEBUG_CHEM_)
                    {
                        fmt::print("Set total concentration zero for {}.\n", buffer_name_[size_rock_]);
                    }
                    c_buffer_[size_rock_] = 0.;
                }
                else
                {
                    c_buffer_[size_rock_] = std::stod(lineTab[0]);

                    if (lineTab.size() == 1) {
                        if (PRINT_DEBUG_CHEM_)
                        {
                            fmt::print("Set Log activity to zero for {}.\n", buffer_name_[size_rock_]);
                        }
                        log_af_[size_rock_] = 0.;
                    }
                    else
                    {
                        log_af_[size_rock_] = std::stod(lineTab[1]);
                        if (lineTab.size() == 2)
                        {
                            if (PRINT_DEBUG_CHEM_)
                            {
                                fmt::print("{}  is in aqueous phase\n", buffer_name_[size_rock_]);
                            }
                        }
                        else
                        {
                            std::string NAME = to_upper_case(lineTab[2]);
                            if (NAME == "OIL")
                            {
                                phase_[size_rock_] = gFluidPhase::OIL;
                                pos_gas_[size_gas_] = size_rock_;
                                ++size_gas_;
                            }
                            else if (NAME == "GAS")
                            {
                                phase_[size_rock_] = gFluidPhase::GAS;
                                pos_gas_[size_gas_] = size_rock_;
                                ++size_gas_;
                            }
                        }
                    }
                }
                ++size_rock_;
            }
        }
        return size_rock_ - i0;
    }
}

int InitChem::read_gas_phase(const GeoChemBlockContent& block_content)
{
    if (block_content.empty())
        return 0;
    else
    {
        for (const auto& [buf_name, buf_values] : block_content)
        {
            if (to_upper_case(buf_name) == "VOLUME_FRAC")
            {
                const std::vector<std::string> lineTab = split_string(buf_values);
                if (lineTab.empty())
                {
                    if (PRINT_DEBUG_CHEM_)
                    {
                        fmt::print("Set volume for gas phase equal to water volume\n");
                    }

                    gas_phase_frac_ = 1.;
                }
                else
                {
                    gas_phase_frac_ = std::stod(lineTab[0]);

                }
            } else if (in_list(buf_name, gas_phase_name_) == -1) // not added before
            {
                const std::vector<std::string> lineTab = split_string(buf_values);

                gas_phase_name_.push_back(buf_name);
                if (lineTab.empty())
                {
                    if (PRINT_DEBUG_CHEM_)
                    {
                        fmt::print("Set zero gas pressure for {}.\n", gas_phase_name_[size_gas_phase_]);
                    }
                    gas_phase_pressure_[size_gas_phase_] = 0.;
                }
                else
                {
                    gas_phase_pressure_[size_gas_phase_] = std::stod(lineTab[0]);
                }

                ++size_gas_phase_;
            }
        }
        return size_gas_phase_;
    }
}

int InitChem::read_io(const GeoChemBlockContent& block_content)
{
    if (block_content.empty())
        return 0;
    else
    {
        int i0 = size_io_;
        for (const auto& block_pair : block_content)
        {
            const std::vector<std::string> lineTab = split_string(block_pair.second);

            io_name_.push_back(block_pair.first);
            c_io_[size_io_] = 1e-3;
            if (!lineTab.empty())
                c_io_[size_io_] = std::stod(lineTab[0]);
            size_io_++;
        }
        return size_io_ - i0;
    }
}

int InitChem::read_surface(const GeoChemBlockContent& block_content)
{
    // The order of the keywords will be alphabetical, b/c GeoChemBlockContent
    // is an std::map with strings as keys, and upper-case letters come before
    // lower-case letters. As a consequence, the keyword "METHOD" will be read
    // before "S_AREA", in which case we can end up setting an integration
    // method for the diffuse layer even though we later find out that
    // CALC_DIFFUSE_LAYER == false.
    //
    // This is ok though, since the integration is only triggered when the flag
    // is true.
    //
    if (block_content.empty())
        return 0;

    size_surf_ = 0;
    std::vector<double> pos_surf;
    std::vector<double> c_surf;

    int i0 = size_surf_;
    for (auto it = block_content.begin(); it != block_content.end(); ++it)
    {
        if (in_list((*it).first, surf_name_) == -1)
        {
            const std::vector<std::string> lineTab = split_string((*it).second);

            std::string w = to_upper_case((*it).first);
            if (w == "METHOD")
            {
                if (!lineTab.empty()) surface_options_from_input_.smethod_ = std::stoi(lineTab[0]);
            }
            else if (w == "S_AREA")
            {
                if (!lineTab.empty())
                {
                    SA_ = std::stod(lineTab[0]);
                }
                if (lineTab.size() > 1)
                {
                    d_DL_ = std::stod(lineTab[1]);
                    surface_options_from_input_.CALCULATE_DIFFUSE_LAYER = true;
                }
            }
            else if (w == "VALENCE_SYMMETRICAL"){
                if (!lineTab.empty()){
                    surface_options_from_input_.VALENCE_OF_SYMMETRICAL_ELECTROLYTE = std::stod(lineTab[0]);
                }
            }
            else if (w == "ONLY_COUNTER_IONS")
            {
                surface_options_from_input_.ONLY_COUNTER_IONS = true;
            }
            else if (w == "EXCHANGE")
            {
                surface_options_from_input_.EXCHANGE_MODEL_DL = true;
                GeoChemBlockContent exchange_inp = { {"XDL-", "1e-8"}, {"YDL+", "1e-8"} };
                read_io(exchange_inp);
            }
                // 7/3-22: Replace "Donnan" flag with new keyword "INTEGRATION" that has several options.
                //         If this keyword is not entered, the default method is "AKSEL_HIORTH".
            else if (w == "INTEGRATION")
            {
                if (lineTab.empty())
                {
                    throw ModelDoesNotExistException("Need to provide an integration method for the diffuse layer!");
                }
                else
                {
                    const std::string method = to_upper_case(lineTab[0]);

                    if (method == "DONNAN" || method == "MEAN_POTENTIAL" || method == "MEAN_POTENTIAL_APPROXIMATION")
                    {
                        throw ModelDoesNotExistException("Mean potential approximation is currently not supported!");
                        surface_options_from_input_.INTEGRATION_METHOD = DiffuseLayerIntegrationMethod::MEAN_POTENTIAL_APPROXIMATION;
                    }
                    else if (method == "AH" || method == "AKSEL" || method == "AKSEL_HIORTH")
                    {
                        surface_options_from_input_.INTEGRATION_METHOD = DiffuseLayerIntegrationMethod::AKSEL_HIORTH;
                    }
                    else if (method == "ANALYTICAL" || method == "ANA")
                    {
                        surface_options_from_input_.INTEGRATION_METHOD = DiffuseLayerIntegrationMethod::ANALYTICAL;
                    }
                    else if (method == "BW" || method == "BORKOVEC_WESTALL" || method == "BORKOVEC" || method == "WESTALL")
                    {
                        surface_options_from_input_.INTEGRATION_METHOD = DiffuseLayerIntegrationMethod::BORKOVEC_WESTALL;
                    }
                    else
                    {
                        throw ModelDoesNotExistException("Need to provide an integration method for the diffuse layer!");
                    }
                }
            }
            // Otherwise, assume surface basis specie.
            // Need both name & site density for charged surface:
            else if (!lineTab.empty())
            {
                surf_name_.push_back((*it).first);
                c_surf.push_back(std::stod(lineTab[0]));
                ++size_surf_;
            }
        }
    }

    for (int i = i0; i < size_surf_; ++i)               c_surf_[i] = c_surf[i - i0];

    if(surface_options_from_input_.EXCHANGE_MODEL_DL)   surface_options_from_input_.CALCULATE_DIFFUSE_LAYER = false;

    set_diffuse_layer_options_from_read_input();  // uses value of smethod_

    return size_surf_;
}

bool InitChem::basis_specie_defined(std::size_t pos, const std::vector<std::string>& names) const
{
    const std::string& name      = Basis_db_->row_name_[pos];
    const std::string& nick_name = Basis_db_->nick_name_[pos];

    if(element_is_in_container(name, names))                      return true;
    if(element_is_in_container(nick_name, names))                 return true;
    if(element_is_in_container(to_upper_case(name), names))       return true;
    if(element_is_in_container(to_upper_case(nick_name), names))  return true;
    return false;
}

// Test if the charge of all basis species adds up to the charge of the complex or mineral (zero)
void InitChem::check_db(ChemTable* Test, ChemTable* Basis)
{
    std::vector<double> charge;
    std::vector<double> scharge;
    std::vector<double> charge_tot;
    std::vector<double> basis_exch;

    charge.resize(Test->noRows_);
    scharge.resize(Test->noRows_);
    charge_tot.resize(Test->noRows_);
    basis_exch.resize(Test->noRows_);

    bool all_species_zero   = true;
    bool both_charges       = false;
    bool except             = false;

    for (int i = 0; i < Test->noRows_; ++i)
    {
        charge[i]        = 0.0;
        charge_tot[i]    = 0.0;
        scharge[i]       = 0.0;
        basis_exch[i]    = 0.0;
        all_species_zero = true;

        for (int j = 0; j < Test->noColumns_; ++j)
        {
            double mij = Test->M_[i][j];
            //
            charge_tot[i] += (Basis->charge_[j]+ Basis->scharge_[j]) * mij;
            charge[i]     += Basis->charge_[j]* mij;
            scharge[i]    += Basis->scharge_[j] * mij;
            //
            if(Test->type_[i]==GeochemicalComponentType::ION_EXCHANGE
               && Basis->type_[j]==GeochemicalComponentType::AQUEOUS_COMPLEX)
            {
                basis_exch[i] += Basis->charge_[j] * mij;
            }
            if (std::fabs(Basis->charge_[j]) > 0
                && fabs(Basis->scharge_[j]) > 0)
            {
                both_charges = true;
            }
            if (std::fabs(mij) > 0)
            {
                all_species_zero = false;
            }
        }

        const double dz = charge_tot[i] - (Test->charge_[i]+Test->scharge_[i]);
        if (std::fabs(dz) >= 1e-8)
        {
            if (Test->type_[i] != GeochemicalComponentType::ION_EXCHANGE)
            {
                fmt::print("Check Database, the charge of element {}, does not match the basis specie charge.\n", Test->row_name_[i]);
                except = true;
            }
            else if(Test->charge_[i] != basis_exch[i])
            {   // Ion exchange complex has the charge of basis specie - Gaines Thompson
                fmt::print("Checking database... The charge of element {} does not match the basis specie charge.\n", Test->row_name_[i]);
                except = true;
            }
        }
        if (all_species_zero)
        {
            fmt::print("Checking database... The element {} does not contain any basis species.\n", Test->row_name_[i]);
            except = true;
        }
        if (both_charges)
        {
            fmt::print("Checking database... The element {} contains elements with both surface charge and charge.\n", Test->row_name_[i]);
            except = true;
        }
        if (except)
        {
            fmt::print("{} (charge={}, scharge={}) = ", Test->row_name_[i], Test->charge_[i], Test->scharge_[i]);
            for (std::size_t j = 0; j < Test->sparseM_[i].size(); ++j)
            {
                int x = Test->sparseM_[i][j].first;
                fmt::print("{} {} (charge {}, scharge {})", Test->sparseM_[i][j].second, Basis->row_name_[x], Basis->charge_[x], Basis->scharge_[x]);
            }
            fmt::print("\nEstimated charge=({})\n", charge[i]);
        }
    }
}

void InitChem::set_diffuse_layer_options_from_read_input()
{
    const int surf_method = surface_options_from_input_.smethod_;

    // Note: surface_flag_ also depends on smethod_, but can vary for different
    //       ChemBasVec objects
    //       (depending on whether surface complexes are present or not)
    if(surf_method == 0 || surf_method == 1)
    {
        surface_options_from_input_.GRAHAME_EQ_VERSION = GrahameEquation::ASSUME_CHARGE_BALANCE;
    }
    else if(surf_method == 2)
    {
        surface_options_from_input_.GRAHAME_EQ_VERSION = GrahameEquation::GENERAL;
    }
    else if(surf_method == 3)
    {
        surface_options_from_input_.GRAHAME_EQ_VERSION = GrahameEquation::SYMMETRICAL_ELECTROLYTE;
    }
}

// checks if the map of old basis to new contains duplicates, i.e. SO4-2 and SO4
void InitChem::clean_transformation()
{
    std::vector<std::string> keys, vals;
    for (auto key : basis_transformation_)
    {
        keys.push_back(key.first);
        vals.push_back(key.second);
    }
    for (std::size_t itr = 0; itr < keys.size(); ++itr)
    {
        int bn = Basis_db_->get_row_index(keys[itr]);
        int an = Aq_db_->get_row_index(vals[itr]);
        basis_transformation_.erase(keys[itr]);
        if (bn > -1 && an > -1)
        {
            basis_transformation_[Basis_db_->row_name_[bn]] = Aq_db_->row_name_[an];
        } // if not found we remove the basis transformation
        else {
            if(PRINT_DEBUG_CHEM_)
            fmt::print("Warning: Basis specie {} cannot be exchanged by {}", keys[itr], vals[itr]);
        }
    }



    for (auto key : basis_transformation_)
    {
        int bn = Basis_db_->get_row_index(key.first);

        for (std::size_t i = 0; i < basis_species_name_.size(); ++i)
        {
            int oldb = Basis_db_->get_row_index(basis_species_name_[i]);
            if (oldb == bn)
            {
                basis_species_name_[i] = key.second;
            }
        }
    }

}
void InitChem::basis_transformation()
{
    double** BasisTrans;
    double** BasisTransInv;
    std::vector<int> oldb, newb;
    int sizeB = Basis_db_->noRows_;
    BasisTrans = allocateMemoryForMatrix(sizeB);
    BasisTransInv = allocateMemoryForMatrix(sizeB);

    std::vector<double> mol_weight;
    for (auto key : basis_transformation_)
    {
        int bn = Basis_db_->get_row_index(key.first);
        int an = Aq_db_->get_row_index(key.second);
        oldb.push_back(bn);
        newb.push_back(an);
    }
    assert(sizeB == (Aq_db_->noColumns_-ChemTable::skip_col_complex()));
 //   for (auto i : sort_vector_indices(oldb,/*reverse=*/ true))
 //   {

 //   }
    for (int i = 0; i < sizeB; ++i)
        for (int j = 0; j < sizeB; ++j)
            BasisTrans[i][j] = (i==j?1.:0.);

    for (std::size_t i = 0; i < oldb.size(); ++i)
    {
        double Mw = 0.;
        for (int j = 0; j < sizeB; ++j)
        {
            BasisTrans[oldb[i]][j] = Aq_db_->M_[newb[i]][ChemTable::skip_col_complex() + j];
            Mw += Aq_db_->M_[newb[i]][ChemTable::skip_col_complex() + j] * Basis_db_->mol_weight_[j];
        }
        mol_weight.push_back(Mw);
    }
    double det;
    invert_matrix(BasisTrans, sizeB, BasisTransInv, det);


    // Basis_db_   : Name	NickName	type	charge	scharge	mol_weight	a0	DeltaG	DeltaH	S	a1	a2	a3	a4	c1	c2	omega
    // Aq_db_      : Name	NickName	type	charge	scharge	            a0	DeltaG	DeltaH	S	a1	a2	a3	a4	c1	c2	omega
    // Mineral_db_ : Name	NickName	mol_volume	DeltaG	DeltaH	S	a1	a2	a3	a7

    std::map<std::string, std::vector<double >&>::iterator it;
    std::map<std::string, std::vector<std::string >&>::iterator its;
    for (std::size_t i = 0; i < oldb.size(); ++i)
    {
        for (it = Aq_db_->double_vector_map_.begin(); it != Aq_db_->double_vector_map_.end();++it)
        {
            double value_old,value_new;
            bool found_old = Basis_db_->get_value(it->first, oldb[i], value_old);
            bool found_new = Aq_db_->get_value(it->first, newb[i], value_new);
            if (found_old && found_new) //swap
            {
                Basis_db_->set_value(it->first, oldb[i], value_new);
                Aq_db_->set_value(it->first, newb[i], value_old);
            }
            if (it->first == "mol_weight") // mol weight only for basis species
            {
                Basis_db_->set_value(it->first, oldb[i], mol_weight[i]);
            }
        }

        for (its = Aq_db_->string_vector_map_.begin(); its != Aq_db_->string_vector_map_.end(); ++its)
        {
            std::string value_old, value_new;
            bool found_old = Basis_db_->get_value(its->first, oldb[i], value_old);
            bool found_new = Aq_db_->get_value(its->first, newb[i], value_new);
            if (found_old && found_new) //swap
            {
                Basis_db_->set_value(its->first, oldb[i], value_new);
                Aq_db_->set_value(its->first, newb[i], value_old);
            }
        }

    }
    // add correct matrix elements for aq complex, e.g. Ca+2 prev basis element now aqueous complex
    for (std::size_t i = 0; i < oldb.size(); ++i)
    {
        for (int j = 0; j < sizeB; ++j)
        {
            Aq_db_->M_[newb[i]][ChemTable::skip_col_complex()+j] = 0.;
        }
        Aq_db_->M_[newb[i]][ChemTable::skip_col_complex()+oldb[i]] = 1.;
    }
    // Here we do the transformation
    std::vector<double> m_elem(sizeB); // TODO: use sparse matrices
    for (int row = 0; row < Aq_db_->noRows_; ++row)
    {
        for (int i = 0; i < sizeB; ++i)
        {
            m_elem[i] = 0.;
            for (int j = 0; j < sizeB; ++j)
            {
                m_elem[i] += Aq_db_->M_[row][ChemTable::skip_col_complex() + j]* BasisTransInv[j][i];
            }
        }
        for (int i = 0; i < sizeB; ++i) Aq_db_->M_[row][ChemTable::skip_col_complex()+i] = m_elem[i];
    }
    for (int row = 0; row < Mineral_db_->noRows_; ++row)
    {
        for (int i = 0; i < sizeB; ++i)
        {
            m_elem[i] = 0.;
            for (int j = 0; j < sizeB; ++j)
            {
                m_elem[i] +=  Mineral_db_->M_[row][ChemTable::skip_col_min() + j]* BasisTransInv[j][i];
            }
        }
        for (int i = 0; i < sizeB; ++i) Mineral_db_->M_[row][ChemTable::skip_col_min() + i] = m_elem[i];
    }

    for (auto key : basis_transformation_) // change name of basis species
    {
        int col=Aq_db_->get_col_index(key.first);
        Aq_db_->col_name_[col] = key.second;
        col = Mineral_db_->get_col_index(key.first);
        Mineral_db_->col_name_[col] = key.second;
    }


    if (BasisTrans)        freeMatrixMemory(BasisTrans, sizeB);
    if (BasisTransInv)   freeMatrixMemory(BasisTransInv, sizeB);
}

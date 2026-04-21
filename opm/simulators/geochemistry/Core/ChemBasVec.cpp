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
#include <opm/simulators/geochemistry/Core/ChemBasVec.h>
#include <opm/simulators/geochemistry/Core/ChemGCSolver.h>
#include <opm/simulators/geochemistry/Core/ChemInitChem.h>

BasVec* BasVec::createFromInitChem(InitChem* ICS)
{
    return new BasVec(ICS);
}

BasVec::BasVec(InitChem* ICS_in)
    : ICS_(ICS_in)
    , charge_balance_(ICS_in->charge_balance_)
    , surface_options_(ICS_in->surface_options_from_input_)
    , WHTO_(ICS_in->WHTO_)
    , d_DL_(ICS_in->d_DL_)
    , SA_(ICS_in->SA_)
{
    std::fill(Temp_.begin(), Temp_.end(), ICS_->Temp_);
    std::fill(Pres_.begin(), Pres_.end(), ICS_->Pres_);

    // Subset of full databases
    ChemTable& SM_basis = *ICS_->SM_basis_;
    ChemTable& SM_all = *ICS_->SM_all_;
    ChemTable& SM_mineral = *ICS_->SM_mineral_;

    const int full_basis_size = SM_all.noColumns_;
    const int full_buffer_size = SM_mineral.noRows_;

    // Allocate memory & initialize with sensible defaults
    size_ = full_basis_size;

    ctot_.resize(size_, 0.0);
    ctot_ads_.resize(size_, 0.0);
    ctot_dl_excess_.resize(size_, 0.0);
    delta_mineral_.resize(size_, 0.0);

    ctot_mineral_.resize(full_buffer_size, 0.0);

    log_m_.resize(size_, NumericalConstants::LOG10_ALMOST_ZERO);
    log_a_.resize(size_, NumericalConstants::LOG10_ALMOST_ZERO);
    log_g_.resize(size_, 0.0);
    log_a_gas_.resize(ICS_->size_gas_, NumericalConstants::LOG10_ALMOST_ZERO);

    phase_.resize(size_, -1);

    pos_mass_.resize(size_, 0);
    pos_buffer_.resize(size_, 0);
    pos_rock_.resize(size_, 0);
    pos_min_.resize(size_, -1);
    pos_sup_min_.resize(size_, -1);
    pos_min_bas_.resize(size_, -1);
    pos_sup_bas_.resize(size_, -1);
    pos_gas_.resize(size_, -1);
    pos_gas_phase_.resize(ICS_->size_gas_phase_, -1);
    pos_rel_sup_min_.resize(size_, -1);
    pos_rel_gas_.resize(size_, -1);

    gas_phase_pressure_.resize(ICS_->size_gas_phase_, -1);
    c_gas_phase_.resize(ICS_->size_gas_phase_, -1);
    mol_fraction_gas_phase_.resize(ICS_->size_gas_phase_, -1);
    rate_laws_.resize(ICS_->size_sup_min_);

    log_bas_.resize(size_, 0.0);

    for (auto& ctot_vec: ctot_calc_)  ctot_vec.resize(size_, 0.0);

    // Special species H20, E, H+
    pos_water_ = SM_basis.get_row_index("H2O");
    pos_exp_ = SM_basis.get_row_index("E");
    pos_pH_ = SM_basis.get_row_index("H");
    pos_pe_ = SM_basis.get_row_index("E-");

    bool basisIncludesSurfaceComplexes = false;

    for(int basisIdx=0; basisIdx < size_; ++basisIdx)
    {
        if (SM_basis.type_[basisIdx] == GeochemicalComponentType::ION_EXCHANGE) has_ion_exchange_ = true;

        // What is the purpose of the last condition? Can we have E present without any actual surface complexes??
        if(SM_basis.type_[basisIdx] == GeochemicalComponentType::SURFACE_COMPLEX && pos_exp_ != basisIdx)
        {
            basisIncludesSurfaceComplexes = true;
        }
    }

    // The handling of model options such as this could still be improved...
    if(basisIncludesSurfaceComplexes)
    {
        // Need to decide if an extra equation is to be included in the non-linear system
        const int s_method = ICS_->surface_options_from_input_.smethod_;

        if(s_method == 0 || s_method > 3) surface_flag_ = 0;
        else surface_flag_  = 1;  // methods 1, 2, 3
    }
    else
    {
        // Surface chemistry options have previously been set from ICS_ in the constructor initializer list.
        // However, there are not any surface complexes for this particular basis, so we set
        // the diffuse layer flag to false regardless of the value stored in InitChem...
        surface_options_.CALCULATE_DIFFUSE_LAYER = false;
    }

    for (int i = 0; i < ICS_->size_basis_; ++i) // need to set c_tot_ before assigning default values
    {
        const int pos = ICS_->pos_[i];
        ctot_[pos] = ICS_->c_vchem_[i];
    }
    // Initialize activities and molarities
    assert(pos_water_ > -1); // Must always have H2O.
    log_a_[pos_water_] = 0.0;
    log_m_[pos_water_] = 0.0;

    assert(pos_pH_ > -1);  // Must always have H+
    if (ctot_[pos_pH_] > WHTO_)
    {
        log_a_[pos_pH_] = std::log10(ctot_[pos_pH_] - WHTO_);
        log_m_[pos_pH_] = log_a_[pos_pH_];
    }
    else{
        const double est_log_a = std::log10(ctot_[pos_pH_]);
        log_a_[pos_pH_] = (-est_log_a<1.0e-8 || -est_log_a >= 20.0 ? -7.:est_log_a);
        log_m_[pos_pH_] = log_a_[pos_pH_];
    }

    if(pos_exp_ > -1)  // In theory, E could appear at index 0
    {
        log_a_[pos_exp_] = 0.0;
        log_m_[pos_exp_] = 0.0;
    }
    if (pos_pe_ > -1)
    {
        log_a_[pos_pe_] = std::log10(ctot_[pos_pe_]);
        log_m_[pos_pe_] = std::log10(ctot_[pos_pe_]);
    }
    for(int i=0; i < ICS_->size_basis_; ++i)
    {
        const int pos = ICS_->pos_[i];
        if(pos == pos_water_ || pos == pos_pH_ || pos == pos_exp_)
        {
            continue;  // Already initialized
        }

        ctot_[pos] = std::max(1.0e-8, ICS_->c_vchem_[i]);
        log_a_[pos] = std::log10(ICS_->c_vchem_[i]);
        log_m_[pos] = log_a_[pos];
        ctot_ads_[pos] = ICS_->c_ads_[i];
    }

    // Ion exchange species
    for (int i = 0; i < size_; ++i)
    {
        if (SM_basis.type_[i] == GeochemicalComponentType::ION_EXCHANGE)
        {
            for (int j = 0; j < SM_all.noRows_; ++j)
            {
                if (SM_all.M_[j][i] != 0)
                {
                    // Note: We store multiple copies of the same basis specie index (e.g., for X-).
                    pos_exch_bas_.push_back(i);
                    pos_exch_comp_.push_back(j);
                }
            }
        }
    }

    assign_basis2mineral();

    size_gas_ = 0;
    size_mass_ = 0;
    size_rock_ = 0;

    // Buffers
    int count_gas = 0;
    int count_buffer = 0;

    for(int j=0; j < ICS_->size_rock_; ++j)
    {
        for(int i=0; i < SM_mineral.noRows_; ++i)
        {
            if(SM_mineral.name_is_in_row(ICS_->buffer_name_[j], i))
            {
                if(ICS_->c_buffer_[j] > 0.0)
                {
                    const auto phase_type = ICS_->phase_[j];
                    phase_[count_buffer] = phase_type;
                    pos_buffer_[count_buffer] = i;
                    if (phase_type == gFluidPhase::OIL || phase_type == gFluidPhase::GAS)
                    {
                        // Oil and gas buffers are treated the same
                        pos_gas_[count_gas]     = i;
                        pos_rel_gas_[count_gas] = count_buffer;
                        log_a_gas_[count_gas]   = ICS_->log_af_[j];
                        ++count_gas;
                    }
                    ++count_buffer;
                }
            }
        }
    }
    size_gas_ = count_gas;

    // Check if all the basis species are valid species
    int count_mass = 0;
    int count_rock = 0;

 //   for(const auto& specie_name: ICS_->basis_species_name_)
    for (auto specie_name : ICS_->basis_species_name_)
    {
        for(int i=0; i < size_; ++i)
        {
            if(SM_basis.name_is_in_row(specie_name, i))
            {
                // Check if basis specie has an associated buffer...
                const int el = index_of_element_in_array(i, ICS_->pos_rock_.data(), ICS_->size_rock_);

                if(el > -1 && ICS_->c_buffer_[el] > 0.0)  // Basis switching
                {
                    ctot_mineral_[pos_buffer_[el]] = ICS_->c_buffer_[el];
                    pos_rock_[count_rock] = i;
                    ++count_rock;
                }
                else if (specie_name == "E-")
                {
                    continue;
                }
                else if (specie_name == "H" && charge_balance_)
                {
                    continue;
                }
                else if (specie_name == "H" && ICS_->keep_ph_fixed_)
                {
                    continue;
                }
                else // Mass balance is used
                {
                        pos_mass_[count_mass] = i;
                        ++count_mass;
                }
            }
        }
    }
    size_mass_ = count_mass;
    size_rock_ = count_rock;

    // Gas phase
    gas_phase_frac_ = ICS_->gas_phase_frac_;
    size_gase_phase_ = ICS_->size_gas_phase_;
    double moles_gas_tot = 0.;
    for (int i = 0; i < size_gase_phase_; ++i)
    {
        pos_gas_phase_[i] = ICS_->pos_gas_phase_[i];
        gas_phase_pressure_[i] = ICS_->gas_phase_pressure_[i];
        // initially we calculate moles or concentration from ideal gas law always (consistent with PHREEQC)
        c_gas_phase_[i] = gas_phase_pressure_[i] * gas_phase_frac_ / PhysicalConstants::IdealGasConstant / Temp_[gFluidPhase::GAS]*1e-3;// mol/L
        moles_gas_tot += c_gas_phase_[i];
    }
    for (int i = 0; i < size_gase_phase_; ++i)
    {
        mol_fraction_gas_phase_[i] = c_gas_phase_[i] / moles_gas_tot;
        pos_all_gas_buffers_.push_back(pos_gas_phase_[i]);
        size_all_gas_buffers_++;
    }
    add_gas_phase_to_total_conc();
    // gas phases present in gas_phase_ treated seperately

    for (int i = 0; i < static_cast<int>(ICS_->SM_mineral_->row_name_.size()); ++i)
    {
        const int el = index_of_element_in_array(i, ICS_->pos_gas_phase_.data(), ICS_->size_gas_phase_);
        if (el < 0)
        {
            if (string_contains(ICS_->SM_mineral_->row_name_[i], "(G)") || string_contains(ICS_->SM_mineral_->row_name_[i], ",g"))
            {
                pos_all_gas_buffers_.push_back(i);
                size_all_gas_buffers_++;
            }
        }
    }
    // For each rock species, there has to be a corresponding buffer.
    if(count_buffer != size_rock_)
    {
        error_and_exit(
            "Non-matching dimensions for rock and buffer basis species: #rock basis species={:d}, #buffer basis species={:d}.",
            size_rock_,
            count_buffer);
    }

    // VERY IMPORTANT: The following code ASSUMES that minerals are linearly independent.
    assert(number_of_non_null_values(ICS_->rate_) == static_cast<std::size_t>(ICS_->size_sup_min_));
    size_sup_min_ = 0;
    for(int i=0; i < ICS_->size_sup_min_; ++i)
    {
        if(ICS_->c_sup_min_[i] > 0)
        {
            pos_rel_sup_min_[size_sup_min_]         = i;
            pos_sup_bas_[size_sup_min_]             = ICS_->pos_sup_bas_[i];
            pos_sup_min_[size_sup_min_]             = ICS_->pos_sup_min_[i];
            ctot_mineral_[ICS_->pos_sup_min_[i]]    = ICS_->c_sup_min_[i];

            assert(ICS_->rate_[i].has_value());
            rate_laws_[size_sup_min_] = ICS_->rate_[i].value();
            ++size_sup_min_;
        }

        ICS_->SM_mineral_->log_af_[ICS_->pos_sup_min_[i]] = ICS_->log_af_sup_[i];
    }

    // Find reduced stoichiometric matrix for super-saturated minerals
    size_min_ = size_sup_min_ + size_rock_;
    if (size_min_ > 0)
    {
        for (int i = 0; i < size_sup_min_; ++i)
        {
            pos_min_[i] = pos_sup_min_[i];
            pos_min_bas_[i] = pos_sup_bas_[i];
        }

        for (int i = 0; i < size_rock_; ++i)
        {
            pos_min_[size_sup_min_ + i] = pos_buffer_[i];
            pos_min_bas_[size_sup_min_ + i] = pos_rock_[i];
        }
        sm_min_ = allocateMemoryForMatrix(size_min_);
        update_sm_buf(size_min_, pos_min_.data(), pos_min_bas_.data(), sm_min_);
    }
    update_basis_buffer_matrices();
}

BasVec::~BasVec()
{
    if (beta_bas_inv_) freeMatrixMemory(beta_bas_inv_, size_);
    if (sm_buf_) freeMatrixMemory(sm_buf_, size_rock_);
    if (sm_min_) freeMatrixMemory(sm_min_, size_min_);

    // More than one Vchem points to a single ICS_, should be freed separately.
    ICS_ = nullptr;

    delete Vp_eq_;
    Vp_eq_ = nullptr;
}

double BasVec::update_exchange_capacity_dl_model(double frac_DL)
{
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    const int posXDL = SM_basis.get_row_index("XDL-");
    const int posYDL = SM_basis.get_row_index("YDL+");

    double CEC = -WHTO_ * frac_DL;  // to get correct H+ concentration
    double AEC = 0.0;
    for (int i = 0; i < size_; ++i)
    {
        if (SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            if (SM_basis.charge_[i] > 0)
            {
                CEC += SM_basis.charge_[i] * ctot_calc_[GeochemicalComponentType::AQUEOUS_COMPLEX][i] * frac_DL;
            }
            else if (SM_basis.charge_[i] < 0)
            {
                AEC += -SM_basis.charge_[i] * ctot_calc_[GeochemicalComponentType::AQUEOUS_COMPLEX][i] * frac_DL;
            }
        }
    }

    if (sch_ < 0.0)     CEC += -sch_ * SA_ / PhysicalConstants::Faraday;
    else                AEC += sch_ * SA_ / PhysicalConstants::Faraday;

    // odn (12/7-22): To avoid potential negative CEC-value due to H+-term.
    ctot_[posXDL] = std::max(NumericalConstants::ZERO_CONCENTRATION_THRESHOLD, CEC);
    ctot_[posYDL] = AEC;

    if (sch_ < 0.0)     return CEC;
    else                return AEC;
}

int BasVec::dimension() const
{
    return size_mass_ + surface_flag_ + size_gas_;
}

/**
 * In addition to updating total concentrations inside the current BasVec
 * instance, store changes in the mineral concentrations in the out-parameter.
 *
 * @param [in] dc Changes in basis species concentrations over the last time
 *                step (old - new). Accounts for both aqueous species
 *                and surface species.
 * @param [out] dc_min Computed changes in the mineral concentrations.
 */
void BasVec::update_ctot_and_ctot_mineral(const std::vector<double>& dc,
                                          std::vector<double>& dc_min)
{
    for (int i = 0; i < ICS_->size_basis_; ++i)
    {
        ctot_[ICS_->pos_[i]] = dc[ICS_->kcmap_[i]];
    }

    for (int i = 0; i < ICS_->size_sup_min_; ++i)
    {
        ctot_mineral_[ICS_->pos_sup_min_[i]] = 0.0;
    }

    for (int i = 0; i < size_min_; ++i)
    {
        double c_tmp = 0.0;
        for (int j = 0; j < size_min_; ++j)
        {
            c_tmp += sm_min_[i][j] * ctot_[pos_min_bas_[j]];
        }
        ctot_mineral_[pos_min_[i]] = c_tmp;
    }

    for (int i = 0; i < ICS_->size_min_; ++i)
    {
        dc_min[ICS_->krmap_[i]] = ctot_mineral_[ICS_->pos_min_[i]];
    }
}

void BasVec::write(const std::string& name)
{
    const ChemTable& SM_mineral = *ICS_->SM_mineral_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    FILE* fpd = my_fopen(name.c_str(), "w");
    fmt::print(fpd, "Temp[K]={:g}\t Pressure [bar] {:g}\n", Temp_[gFluidPhase::WATER], Pres_[gFluidPhase::WATER]);
    fmt::print(fpd,"name\t c_tot\t log_m\t log_a\t log_g\t c_ads\n");
    for(int i=0; i < size_; ++i)
    {
        fmt::print(
            fpd,
            "{:s}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\n",
            SM_basis.row_name_[i],
            ctot_[i],
            log_m_[i],
            log_a_[i],
            log_g_[i],
            ctot_ads_[i]);
    }
    fmt::print(fpd,"Mass balance:\n");
    for(int i=0; i < size_mass_; ++i) fmt::print(fpd,"{:s} ", SM_basis.row_name_[pos_mass_[i]]);
    fmt::print(fpd, "\nRock basis:\n");
    for(int i=0; i < size_rock_; ++i) fmt::print(fpd,"{:s} ", SM_basis.row_name_[pos_rock_[i]]);
    fmt::print(fpd, "\nSup basis:\n");
    for(int i=0; i < size_sup_min_; ++i) fmt::print(fpd, "{:s} ",SM_basis.row_name_[pos_sup_bas_[i]]);
    fmt::print(fpd, "\nRock buffer:\n");
    for(int i=0; i < size_rock_; ++i) fmt::print(fpd, "{:s} ", SM_mineral.row_name_[pos_buffer_[i]]);
    fmt::print(fpd, "\nSupersaturated minerals:\n");
    for(int i=0; i < size_sup_min_; ++i) fmt::print(fpd, "{:s} ", SM_mineral.row_name_[pos_sup_min_[i]]);
    fmt::print(fpd, "\n");

    if (charge_balance_)
    {
        fmt::print(fpd, "pH determined by charge balance\n");
    }

    if (size_rock_ > 0)
    {
        assert(beta_bas_inv_ != nullptr);

        fmt::print(fpd, "Equilibrium: Basis transformation matrix:\n\t");
        for (int i = 0; i < size_; ++i)
        {
            int posi = index_of_element_in_array(i, pos_rock_.data(), size_rock_);
            if (posi > -1)
                fmt::print(fpd, "{:s}\t", SM_mineral.row_name_[pos_buffer_[posi]]);
            else
                fmt::print(fpd, "{:s}\t", SM_basis.row_name_[i]);
        }
        fmt::print(fpd, "\n");
        for (int i = 0; i < size_; ++i)
        {
            fmt::print(fpd, "{:s}\t", SM_basis.row_name_[i].c_str());
            for (int j = 0; j < size_; ++j) fmt::print(fpd, "{:g}\t", beta_bas_inv_[i][j]);
            fmt::print(fpd, "\n ");
        }
    }

    if (size_min_ > 0)
    {
        fmt::print(fpd, "Relation between basis fluxes and mineral fluxes:\n");
        for (int i = 0; i < size_min_; ++i)
        {
            fmt::print(fpd, "{:s} = ", SM_mineral.row_name_[pos_min_[i]].c_str());
            int jj = -1;
            for (int j = 0; j < size_min_; ++j)

            {
                double ff = fabs(sm_min_[i][j]);
                if (ff > 0)
                {
                    ++jj;
                    if (sm_min_[i][j] > 0 && jj > 0) fmt::print(fpd, "+ {:g}", ff);
                    else if (sm_min_[i][j] < 0 && jj>0) fmt::print(fpd, "- {:g}", ff);
                    else fmt::print(fpd, " {:g}", ff);
                    fmt::print(fpd, "{:s} ", SM_basis.row_name_[pos_min_bas_[j]].c_str());
                }
            }
            fmt::print(fpd, "\n ");
        }
    }
    fclose(fpd);
}

// write solution chemistry in tabulated form
void BasVec::write_solution_chemistry(const std::string& name,bool include_reactions)
{
    [[maybe_unused]] const ChemTable& SM_mineral = *ICS_->SM_mineral_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;
    const ChemTable& SM_all = *ICS_->SM_all_;
    std::vector<double> log_m, log_a, log_g,logK,mol_volume;
    std::vector<std::string> names= SM_all.row_name_;
    std::vector<std::string> reactions;
    log_m = SM_all.log_m_; log_m.insert(log_m.end(), log_m_.begin(), log_m_.end());
    log_a = SM_all.log_a_; log_a.insert(log_a.end(), log_a_.begin(), log_a_.end());
    log_g = SM_all.log_g_; log_g.insert(log_g.end(), log_g_.begin(), log_g_.end());
    logK = SM_all.logK_; logK.insert(logK.end(), SM_basis.logK_.begin(), SM_basis.logK_.end());
    mol_volume = SM_all.mol_volume_; mol_volume.insert(mol_volume.end(), SM_basis.mol_volume_.begin(), SM_basis.mol_volume_.end());
    names = SM_all.row_name_; names.insert(names.end(), SM_basis.row_name_.begin(), SM_basis.row_name_.end());
    if (include_reactions)
    {
        for (int i = 0; i < SM_all.noRows_; ++i)
            reactions.push_back(SM_all.get_reaction(i));
        for (int i = 0; i < SM_basis.noRows_; ++i)
            reactions.push_back("-");
    }
    FILE* fpd = my_fopen(name.c_str(), "w");
    fmt::print(fpd, "name\tlog_m\tlog_a\tlog_g\tlogK\tmol_volume");
    if (include_reactions)
        fmt::print(fpd, "\treaction");
    fmt::print(fpd, "\n");
    for (int i : sort_vector_indices(log_m))
    {
        fmt::print(
            fpd,
            "{:s}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}",
            names[i],
            log_m[i],
            log_a[i],
            log_g[i],
            logK[i],
            mol_volume[i]);
        if (include_reactions)
            fmt::print(fpd, "\t{:s}", reactions[i]);
        fmt::print(fpd, "\n");
    }

    fclose(fpd);
}

void BasVec::write_buffers(const std::string& name,  bool include_reactions)
{
    FILE* fpd = my_fopen(name.c_str(), "w");
    const ChemTable& SM_mineral = *ICS_->SM_mineral_;
    std::vector<std::string> reactions;
    if (include_reactions)
    {
        for (int i = 0; i < SM_mineral.noRows_; ++i)
            reactions.push_back(SM_mineral.get_reaction(i));
    }
    fmt::print(fpd, "name\tSI\tlogK");
    if (include_reactions)
        fmt::print(fpd, "\treaction");
    fmt::print(fpd, "\n");

    for (int i : sort_vector_indices(SM_mineral.log_a_))
    {
        fmt::print(fpd, "{:s}\t{:g}\t{:g}", SM_mineral.row_name_[i], SM_mineral.log_a_[i], SM_mineral.logK_[i]);
        if (include_reactions)
            fmt::print(fpd, "\t{:s}", reactions[i]);
        fmt::print(fpd, "\n");
    }
    fclose(fpd);
}
void BasVec::write_gas_phase_solution(const std::string& name)
{
    FILE* fpd = my_fopen(name.c_str(), "w");
    const ChemTable& SM_mineral = *ICS_->SM_mineral_;
    fmt::print(fpd,"name\tmol_fraction\tpartial_pressure\tfugacity\n");
    for (int i = 0; i < size_gase_phase_; ++i)
    {
        fmt::print(fpd, "{:s}\t{:g}\t{:g}\t{:g}", SM_mineral.row_name_[pos_gas_phase_[i]], mol_fraction_gas_phase_[i], gas_phase_pressure_[i], SM_mineral.fugacity_[pos_gas_phase_[i]]);
        fmt::print(fpd, "\n");
    }
}
double BasVec::get_species_concentration(const std::string& species_name) const
{
    const ChemTable& SM_basis = *ICS_->SM_basis_;
    const ChemTable& SM_all = *ICS_->SM_all_;

    const int basisIndex = SM_basis.get_row_index(species_name);
    if(basisIndex > -1){
        return std::pow(10.0, log_m_[basisIndex]);
    }

    // If the current species is not a basis species, continue the search...
    const int secondaryIndex = SM_all.get_row_index(species_name);
    if(secondaryIndex > -1){
        return std::pow(10.0, SM_all.log_m_[secondaryIndex]);

    }
    return 0.0;  // If non-existent species, set concentration to zero
}

double BasVec::get_delta_mineral(const std::string& mineral_name) const
{
    if(size_rock_ == 0)     return 0.0;

    const ChemTable& SM_mineral = *ICS_->SM_mineral_;
    const int bufferIndex = SM_mineral.get_row_index(mineral_name);
    if(bufferIndex > -1){
        return SM_mineral.delta_[bufferIndex];
    }
    return 0.0;
}

void BasVec::update_before_write_solution()
{
    if(ICS_->SM_mineral_) calc_complex_sparse(*ICS_->SM_mineral_);
    calc_ctot_aq();
}

// call after write_solution
void BasVec::write_tot_conc(const std::string& name)
{
    FILE* fp = my_fopen(name.c_str(), "w");
    fmt::print(fp, "Species\tch\tc_tot\tc_tot_aq\tc_tot_SC\tc_tot_IO\tc_dl\n");
    const ChemTable& SM_basis = *ICS_->SM_basis_;
    for (int i = 0; i < size_; ++i)
    {

        fmt::print(fp,
            "{:s}\t{:g}\t{:12.8e}\t{:12.8e}\t{:12.8e}\t{:12.8e}\t{:12.8e}\n",
            SM_basis.row_name_[i].c_str(),
            SM_basis.charge_[i],
            ctot_[i],
            ctot_calc_[0][i],
            ctot_calc_[1][i],
            ctot_calc_[2][i],
            ctot_dl_excess_[i]);

    }
    fclose(fp);

}
void BasVec::write_solution(const std::string& name, DiffuseLayerProperties* EDL_props)
{
    // Do some necessary calculations...
    update_before_write_solution();

    FILE* fp = my_fopen(name.c_str(), "w");

    const ChemTable& SM_all = *ICS_->SM_all_;
    const ChemTable& SM_mineral = *ICS_->SM_mineral_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    double ch_DL = 0.0;
    for (int i = 0; i < size_; ++i) ch_DL += SM_basis.charge_[i] * ctot_dl_excess_[i];

    const auto phaseIdx = gFluidPhase::WATER;
    fmt::print(fp, "pH\t Sch[C/m^2]\tSch_DL[C/m^2]\t psi[mV]\t F[C/mol]\t S[m^2/L]\tcbal[eq/L]\n");
    fmt::print(
        fp,
        "{:g}\t {:g}\t {:g}\t {:g}\t {:g}\t {:g}\t{:12.16e}\n",
        -log_a_[pos_pH_],
        sch_,
        ch_DL* PhysicalConstants::Faraday/SA_,
        1.0e3*psi_,
        PhysicalConstants::Faraday,
        SA_,
        calc_solution_charge()
    );
    fmt::print(fp, "Temp\tPres[bar]\tIo_\tpermittivity\tdensity_water[kg/m^3]\n");
    fmt::print(
        fp,
        "{:g}\t{:g}\t{:g}\t{:g}\t{:g}\n",
        Temp_[phaseIdx] - 273.15,
        1.0e-5*Pres_[phaseIdx],
        Io_,
        ICS_->CP_.ew_,
        ICS_->CP_.rho_w_
    );

    if (EDL_props)
    {
        fmt::print(fp, "n(1)\tn(-1)\tn(2)\tn(-2)\tg(1)\tg(-1)\tg(2)\tg(-2)\tfrac_DL\n");
        fmt::print(fp, "{:g}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\n", EDL_props->n1_, EDL_props->n_1_, EDL_props->n2_, EDL_props->n_2_, EDL_props->fDL1_, EDL_props->fDL_1_, EDL_props->fDL2_, EDL_props->fDL_2_, EDL_props->frac_DL_);
    }

    fmt::print(fp, "\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:s}\t", SM_basis.row_name_[i]);
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:s}\t", SM_all.row_name_[i]);
    }

    fmt::print(fp, "\ntype\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:d}\t", static_cast<int>(SM_basis.type_[i]));
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:d} \t", static_cast<int>(SM_all.type_[i]));
    }
    fmt::print(fp, "\nlog_a\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:12.16e}\t", log_a_[i]);
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:12.16e}\t", SM_all.log_a_[i]);
    }
    fmt::print(fp, "\nlog_m\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:12.16e}\t", log_m_[i]);
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:12.16e}\t", SM_all.log_m_[i]);
    }
    fmt::print(fp, "\nlog_g\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:12.16e}\t", log_g_[i]);
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:g} \t", SM_all.log_g_[i]);
    }
    fmt::print(fp, "\nlogK\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:g}\t", SM_basis.logK_[i]);
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:g} \t", SM_all.logK_[i]);
    }
    fmt::print(fp, "\ncharge\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:g}\t", SM_basis.charge_[i]);
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:g} \t", SM_all.charge_[i]);
    }
    fmt::print(fp, "\nscharge\t");
    for (int i = 0; i < size_; ++i){
        fmt::print(fp, "{:g}\t", SM_basis.scharge_[i]);
    }
    for (int i = 0; i < SM_all.noRows_; ++i){
        fmt::print(fp, "{:g} \t", SM_all.scharge_[i]);
    }
    fmt::print(fp, "\n");

    fmt::print(fp, "Species\t ch\t c_tot\t c_tot_aq\tc_tot_SC\tc_tot_IO\tc_DLl\n");
    std::vector<double> flux;
    for (int i = 0; i < size_; ++i)
    {
        const double dh = (i == pos_pH_ ? ctot_calc_[0][i] - WHTO_ : ctot_calc_[0][i]);

        fmt::print(fp,
                   "{:s}\t {:g}\t {:12.8e}\t {:12.8e}\t {:12.8e}\t {:12.8e}\t {:12.8e}\n",
                   SM_basis.row_name_[i].c_str(),
                   SM_basis.charge_[i],
                   ctot_[i],
                   dh,
                   ctot_calc_[1][i],
                   ctot_calc_[2][i],
                   ctot_dl_excess_[i]);

        flux.push_back(ctot_[i] - ctot_calc_[0][i] - ctot_calc_[1][i] - ctot_calc_[2][i]);
    }

    if (size_rock_ > 0)
    {
        fmt::print(fp, "------------------------------------------------------------\n");
        std::vector<double> dm_eq;
        dm_eq.resize(size_rock_);
        fmt::print(fp, "Mineral\tmmol/kg-rock\twt%\tmol/kgw\n");
        for (int i = 0; i < size_rock_; ++i)
        {
            dm_eq[i] = 0.;
            for (int j = 0; j < size_rock_; ++j)
            {
                dm_eq[i] += sm_buf_[i][j] * flux[pos_rock_[j]];
            }
            fmt::print(fp, "{:s}\t{:g}\t{:g}\t{:g}\n",
                       SM_mineral.row_name_[pos_buffer_[i]],
                       1.0e3 * dm_eq[i],
                       1.0e2 * dm_eq[i] * SM_mineral.mol_weight_[pos_buffer_[i]],
                       SM_mineral.delta_[pos_buffer_[i]]);
        }
        fmt::print(fp, "------------------------------------------------------------\n");
    }

    fmt::print(fp, "Name\tSI\tlogK\treaction\n");
    for (int i : sort_vector_indices(SM_mineral.log_a_))
    {
        fmt::print(fp, "{:s}\t{:g}\t{:g}\t", SM_mineral.row_name_[i], SM_mineral.log_a_[i], SM_mineral.logK_[i]);
        for (int j = 0; j < SM_mineral.noColumns_; ++j)
        {
            if (fabs(SM_mineral.M_[i][j]) > 0) {
                fmt::print(fp, " {:g}{:s} ", SM_mineral.M_[i][j], SM_basis.row_name_[j]);
            }
        }
        fmt::print(fp, "\n");
    }

    if (has_ion_exchange_)
    {
        fmt::print(fp, "\n");
        fmt::print(fp, "\t Concen-\tEquiv-\tEquivalent\tLog\t \n");
        fmt::print(fp, "Species\ttration\talents\tfraction\tgamma\tlogK\n");
        for (int i = 0; i < SM_all.noRows_; ++i)
        {
            if (SM_all.type_[i] == GeochemicalComponentType::ION_EXCHANGE)
            {
                fmt::print(fp, "{:s}\t {:10.6e}\t{:10.6e}\t{:10.6e}\t{:g}\t{:g}",
                           SM_all.row_name_[i],
                           POW10(SM_all.log_m_[i]),
                           POW10(SM_all.log_m_[i])*fabs(SM_all.charge_[i]),
                           POW10(SM_all.log_a_[i] - SM_all.log_g_[i]),
                           SM_all.log_g_[i],
                           SM_all.logK_[i]);
                fmt::print(fp, "\n");
            }
        }
    }

    if (pos_exp_ >= 0)
    {
        fmt::print(fp, "\n");
        fmt::print(fp, "species\tconcentration\tfractions\tlogK\n");
        for (int i = 0; i < SM_basis.noRows_; ++i)
        {
            double total_site_conc = 1.0;
            if (SM_basis.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
            {
                const auto& surface_complex_name = SM_basis.row_name_[i];

                // FIXME: Here, we assume that there can only be three basis
                //        surface complexes other than E...
                //        If a different surface complexation model than the
                //        current default one is to be used, the next lines
                //        of code should be changed!
                if(anyOf(surface_complex_name, "GCa+", "GCO3-", "GSiOH"))
                {
                    total_site_conc = ctot_[i];
                }
                if (surface_complex_name != "E")
                {
                    fmt::print(fp, "{:s}\t {:10.6e}\t{:g}\t{:g}",
                               SM_basis.row_name_[i],
                               POW10(log_m_[i]),
                               POW10(log_m_[i]) / total_site_conc,
                               SM_basis.logK_[i]);
                    fmt::print(fp, "\n");
                }
            }
        }
        for (int i = 0; i < SM_all.noRows_; ++i)
        {
            if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
            {
                // UPDATE (15/7-22):
                //
                // Before, the same normalizing concentration was used for all secondary surface complexes.
                // This is not correct in general, e.g., if the total number of calcium sites in the
                // carbonate SCM differs from the total number of carbonate sites. Then, GCO3- would always
                // be used to compute fractions, because GCO3- comes after GCa+ in the basis list..
                //
                const double total_site_conc = ctot_[indexOfBasisSurfaceComplex(i)];
                fmt::print(fp, "{:s}\t {:10.6e}\t{:g}\t{:g}",
                           SM_all.row_name_[i],
                           POW10(SM_all.log_m_[i]),
                           POW10(SM_all.log_m_[i]) / total_site_conc,
                           SM_all.logK_[i]);
                fmt::print(fp, "\n");
            }
        }
    }
    fmt::print(fp, "\n");

    fclose(fp);
}


void BasVec::write_info(BasVecInfo* info){

    if(!info){
        return;
    }

    update_before_write_solution();

    std::map<std::string, double> key_solution_properties;
    key_solution_properties["pH"] = -log_a_[pos_pH_];
    key_solution_properties["Surface_charge"] = sch_;
    key_solution_properties["Surface_potential[mV]"] = psi_ * 1000;
    key_solution_properties["Faraday_constant"] = PhysicalConstants::Faraday;
    key_solution_properties["Specific_surface_area[m^2/L]"] = SA_;
    key_solution_properties["Temperature[Celsius]"] = Temp_[gFluidPhase::WATER] - 273.15;
    key_solution_properties["Pressure[bar]"] = 1.0e-5*Pres_[gFluidPhase::WATER];
    key_solution_properties["Ionic_strength"] = Io_;
    key_solution_properties["Water_permittivity"] = ICS_->CP_.ew_;
    key_solution_properties["Water_density"] = ICS_->CP_.rho_w_;
    info->key_solution_properties_ = key_solution_properties;

    std::map<std::pair<std::string, std::string>, double> species_properties;

    const ChemTable& SM_all = *ICS_->SM_all_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    for (int i = 0; i< size_; ++i)
    {
        const std::string basis_species = SM_basis.row_name_[i];
        species_properties.insert({std::make_pair(basis_species, "type"), SM_basis.type_[i]});
        species_properties.insert({std::make_pair(basis_species, "log_a"), log_a_[i]});
        species_properties.insert({std::make_pair(basis_species, "log_m"), log_m_[i]});
        species_properties.insert({std::make_pair(basis_species, "log_g"), log_g_[i]});
        species_properties.insert({std::make_pair(basis_species, "logK"), SM_basis.logK_[i]});
        species_properties.insert({std::make_pair(basis_species, "charge"), SM_basis.charge_[i]});
        species_properties.insert({std::make_pair(basis_species, "scharge"), SM_basis.scharge_[i]});
    }

    for (int i = 0; i< SM_all.noRows_; ++i)
    {
        const std::string species = SM_all.row_name_[i];
        species_properties.insert({std::make_pair(species, "type"), SM_all.type_[i]});
        species_properties.insert({std::make_pair(species, "log_a"), SM_all.log_a_[i]});
        species_properties.insert({std::make_pair(species, "log_m"), SM_all.log_m_[i]});
        species_properties.insert({std::make_pair(species, "log_g"), SM_all.log_g_[i]});
        species_properties.insert({std::make_pair(species, "logK"), SM_all.logK_[i]});
        species_properties.insert({std::make_pair(species, "charge"), SM_all.charge_[i]});
        species_properties.insert({std::make_pair(species, "scharge"), SM_all.scharge_[i]});
    }

    for(int i=0; i < size_; ++i)
    {
        const std::string basis_species = SM_basis.row_name_[i];
        const double species_aqueous_conc = (i == pos_pH_ ? ctot_calc_[0][i] - WHTO_ : ctot_calc_[0][i]);
        species_properties.insert({std::make_pair(basis_species, "c_tot"), ctot_[i]});
        species_properties.insert({std::make_pair(basis_species, "c_tot_aq"), species_aqueous_conc});
        species_properties.insert({std::make_pair(basis_species, "c_tot_SC"), ctot_calc_[1][i]});
        species_properties.insert({std::make_pair(basis_species, "c_tot_IO"), ctot_calc_[2][i]});
        species_properties.insert({std::make_pair(basis_species, "c_dl"), ctot_dl_excess_[i]});
    }

    info->species_properties_ = species_properties;

}

/*
* Calculates total basis species concentrations and stores the updated values
* in ctot_calc_. Also updates delta_mineral_.
*/
void BasVec::calc_ctot_aq()
{
    ChemTable& SM_all = *ICS_->SM_all_;
    const ChemTable& SM_mineral = *ICS_->SM_mineral_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    calc_rock_spec_conc();
    calc_complex_sparse(SM_all);

    // c_tot_j = m_j + sum_i^Nx SM_i,j x_i
    for(int j=0; j < size_; ++j)
    {
        for(int i=0; i < GeochemicalComponentType::NUMBER_OF_NON_MINERAL_TYPES_; ++i){
            ctot_calc_[i][j] = 0.0;
        }

        if(SM_basis.type_[j] != GeochemicalComponentType::ION_EXCHANGE){
            ctot_calc_[SM_basis.type_[j]][j] = POW10(log_m_[j]);
        }

        for(int i=0; i < SM_all.noRows_; ++i){
            ctot_calc_[SM_all.type_[i]][j] += SM_all.M_[i][j]*POW10(SM_all.log_m_[i]);
        }
    }

    // Question (24/12-21): Why is this needed, the type of H+ is always an aqueous complex, no??
    const int index_pH = SM_basis.type_[pos_pH_];
    assert(index_pH == GeochemicalComponentType::AQUEOUS_COMPLEX);
    ctot_calc_[index_pH][pos_pH_] += WHTO_;

    // Store amount of precipitation / dissolution
    for(int i=0; i < size_rock_; ++i)
    {
        const int pos = pos_rock_[i];
        const int pos_buffer = pos_buffer_[i];
        delta_mineral_[pos] = SM_mineral.delta_[pos_buffer];
    }
}

/**
 * Calculates concentrations of rock-buffered species from the basis transformation matrix, i.e.:
 *
 *       log_10 a = beta^-1 * ( log_10 a' + log_10 K)
 *
 * Updates the array log_a_ at relevant positions and, as a consequence, log_m_.
 */
void BasVec::calc_rock_spec_conc()
{
    if (size_rock_ == 0) return;

    // Initialize
    // (used for non-buffered basis species in matrix multiplication below)
    for(int i=0; i < size_; ++i)
    {
        log_bas_[i] = log_a_[i] + ICS_->SM_basis_->logK_[i];
    }

    // Basis species -> basis mineral. This overwrites the initial values.
    for(int i=0; i < size_rock_; ++i)
    {
        const int pos   = pos_rock_[i];
        const int pos_b = pos_buffer_[i];

        // log activity of basis mineral
        log_bas_[pos] = ICS_->SM_mineral_->log_af_[pos_b] + ICS_->SM_mineral_->logK_[pos_b];
    }

    // Finally, do the matrix multiplication
    for(int i=0; i < size_rock_; ++i)
    {
        const int pos = pos_rock_[i];

        log_a_[pos] = 0.0;
        for(int j=0; j < size_; ++j)
        {
            log_a_[pos] += beta_bas_inv_[pos][j]*log_bas_[j];
        }

        log_m_[pos] = log_a_[pos] - log_g_[pos];
    }
}

/** Calculates concentration of complexes (including ion exchange species). */
void BasVec::calc_complex_sparse(ChemTable& sm) const
{
    for (int i = 0; i < sm.noRows_; ++i)
    {
        sm.log_a_[i] = 0.0;
        sm.log_m_[i] = 0.0;

        for(const auto& [index, coeff]: sm.sparseM_[i])
        {
            sm.log_a_[i] += coeff * log_a_[index];
        }
        sm.log_a_[i] -= sm.logK_[i];

        sm.log_m_[i] = sm.log_a_[i] - sm.log_g_[i];

        if (sm.type_[i] == GeochemicalComponentType::ION_EXCHANGE)
        {
            sm.log_m_[i] += sm.log_QK_[i];
        }
    }
}

double BasVec::update_ionic_strength()
{
    const ChemTable& SM_all = *ICS_->SM_all_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    Io_ = 0.0;

    // Basis species
    for(int i=0; i < size_; ++i)
    {
        if(SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX && i != pos_pe_) // skip electron
        {
            Io_ += SM_basis.charge_[i]*SM_basis.charge_[i]*POW10(log_m_[i]);
        }
    }

    // Secondary species
    for(int i=0; i < SM_all.noRows_; ++i)
    {
        if(SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            Io_ += SM_all.charge_[i]*SM_all.charge_[i]*POW10(SM_all.log_m_[i]);
        }
    }
    Io_ *= 0.5;
    return Io_;
}

/*
 * Updates Vchem.
 *
 * @param ICSRock position index to new rock species relative to Vchem.
 * @param ICSBuffer index for new buffer.
 */
void BasVec::add_new_buffer_mineral(int ICSrock, int ICSbuffer)
{
    // Remove the basis species from the mass balance
    std::vector<int> pos_mass2(size_mass_ - 1,  -1);  // temporary workspace

    bool found = false;
    int idum = 0;
    for(int i=0; i < size_mass_; ++i)
    {
        if(pos_mass_[i] == ICSrock){
            found = true;
        }
        else
        {
            pos_mass2[idum] = pos_mass_[i];
            ++idum;
        }
    }

    if(found)
    {
        --size_mass_;
        for(int i=0; i < size_mass_; ++i){
            pos_mass_[i] = pos_mass2[i];
        }
    }
    else
    {
        const auto specie_name = ICS_->SM_basis_->row_name_[ICSrock];
        error_and_exit("BasVec::add_new_buffer_mineral(): Something is wrong, basis specie {:s} is not part of the mass balance..", specie_name);
    }

    // Now, we are ready to add the basis species to the rock buffered ones
    pos_rock_[size_rock_]   = ICSrock;
    pos_buffer_[size_rock_] = ICSbuffer;
    ++size_rock_;
    update_basis_buffer_matrices(/* old dimension= */ size_rock_ - 1);
}

/**
 * Transform rate minerals into equilibrium phases.
 *
 * - Uses charge-balance to calculate H+ (regardless of what is done elsewhere).
 * - Equilibrates with any surface species that may be present.
 *   (Ion exchangers and/or surface complexes.)
 *
 * NOTE:
 *   This function always overwrites previously stored mineral concentrations
 *   with the values given as input. Thus, even if ICS_ has incorrect values
 *   before calling the function, it should not matter as long as the correct
 *   values are passed in.
 *
 * @param [in] cmin Mineral concentrations in mol/L.
 * @return Member variable Vp_eq_.
 *         If it does not already exist, this function creates it.
 */
BasVec* BasVec::convert_nlin_rate_equations_to_equilibrium(const double* cmin)
{
    if(!Vp_eq_)
    {
        const bool old_charge_balance_flag = ICS_->charge_balance_;
        ICS_->charge_balance_ = true;

        for (int i = 0; i < ICS_->size_min_; ++i)
        {
            ICS_->c_mineral_[i] = cmin[ICS_->krmap_[i]];

            if (i < ICS_->size_sup_min_)
            {
                ICS_->c_sup_min_[i] = ICS_->c_mineral_[i];
            }
            else
            {
                if ( ((i - ICS_->size_sup_min_) >= 0) || ( (i - ICS_->size_sup_min_) < ICS_->size_rock_ ) )
                {
                    ICS_->c_buffer_[i - ICS_->size_sup_min_] = ICS_->c_mineral_[i];
                }
                else
                {
                    error_and_exit("Invalid index in convert_nlin_rate_equations_to_equilibrium()...");
                }
            }
        }

        Vp_eq_ = BasVec::createFromInitChem(ICS_);

        if (ICS_->includes_ion_exchange() || ICS_->includes_surface_complexes())
        {
            Vp_eq_->equilibrate_ = true;
        }

        // Initialize
        Vp_eq_->Temp_ = Temp_;

        for(int i=0; i < size_; ++i)
        {
            Vp_eq_->ctot_[i] = ctot_[i];
            Vp_eq_->ctot_ads_[i] = ctot_ads_[i];
            Vp_eq_->log_a_[i] = log_a_[i];
            Vp_eq_->log_m_[i] = log_m_[i];
            Vp_eq_->log_g_[i] = log_g_[i];
        }

        if (Vp_eq_->size_gas_ > 0)
        {
            for (int i = 0; i < Vp_eq_->size_gas_; ++i)
            {
                Vp_eq_->pos_rel_gas_[i] = 0;
                Vp_eq_->pos_gas_[i] = 0;
            }

            Vp_eq_->size_gas_ = 0;
            for (int i = 0; i < Vp_eq_->size_rock_; ++i)
            {
                Vp_eq_->phase_[i] = gFluidPhase::WATER;
            }
        }

        // Add supersaturated mineral to rock buffer for fast equilibration
        // calculation
        for(int i=0; i < size_sup_min_; ++i)
        {
            const int index_rock = Vp_eq_->pos_sup_bas_[i];
            const int index_buffer = Vp_eq_->pos_sup_min_[i];
            assert(index_rock != -1);
            assert(index_buffer != -1);

            fmt::print(
                "Add buffer mineral {:s} with basis {:s}\n",
                Vp_eq_->ICS_->SM_mineral_->row_name_[pos_sup_min_[i]],
                Vp_eq_->ICS_->SM_basis_->row_name_[pos_sup_bas_[i]]);

            Vp_eq_->add_new_buffer_mineral(index_rock, index_buffer);
        }

        Vp_eq_->charge_balance_ = true;
        Vp_eq_->size_sup_min_ = 0;
        ICS_->charge_balance_ = old_charge_balance_flag;  // Set back.
    }
    return Vp_eq_;
}

/** @return -1 if no minerals are supersaturated. */
int BasVec::index_of_most_supersaturated_mineral(const std::vector<int>& mineral_list, int size, const std::vector<int>& old_mineral_list, int old_size) const
{
    // Note: Often, most entries of the input vectors are -1.
    //       We use the size parameters to only loop over the first part (>= 0).
    double a = 1.0e-8;
    int pos_b = -1;
    for(int i=0; i < size; ++i) // checks all mineral phases
    {
        const int pos = mineral_list[i];
        const double b = ICS_->SM_mineral_->log_a_[pos]-ICS_->SM_mineral_->log_af_[pos];

        const bool mineral_not_already_added = (index_of_element_in_array(pos, old_mineral_list.data(), old_size) < 0);
        if(b > a && mineral_not_already_added)
        {
            pos_b = i;
            a = b;
        }
    }

    return pos_b;
}

/**
 * @param dt [in] Time step.
 * @param c_min [out] If there are rate-minerals, the updated concentrations.
 *                    Otherwise, no change in the input concentrations.
 */
void BasVec::update_mineral_concentrations_kinetics(double dt, double* c_min)
{
    if(size_sup_min_ == 0)      return;
    assert(c_min);

    // Otherwise, we are good to go...
    const double aH = POW10(log_a_[pos_pH_]);

    for(int i=0; i < size_sup_min_; ++i)
    {
        const int posi = pos_sup_min_[i];

        double SI = 0.0;
        for(int k=0; k < ICS_->SM_mineral_->noColumns_; ++k)
        {
            // Are minerals always placed first?
            SI += ICS_->SM_mineral_->M_[posi][k]*log_a_[k];
        }
        SI -= ICS_->SM_mineral_->logK_[posi];

        ICS_->SM_mineral_->log_a_[posi] = SI;
        ICS_->SM_mineral_->log_m_[posi] = SI - ICS_->SM_mineral_->log_g_[posi];

        const double sgn = (SI < 0) ? 1.0 : -1.0;
        const double k1_T = rate_laws_[i].k1_;
        const double k2_T = rate_laws_[i].k2_;
        const double m = rate_laws_[i].m_;
        const double n = rate_laws_[i].n_;
        const double n_acid = rate_laws_[i].n_acid_;

        // Update (4/7-22): Rate equation should be consistent with ChemGCSolver::calc_F_sup_min().
        //                  Before, there was no exponent in the pH-term...
        const double fac1 = k1_T + k2_T*pow(aH, n_acid);
        const double fac2 = pow(fabs(1.0-POW10(SI*m)), n);  // POW10(SI) == log10(IAP/K) = Omega.

        ctot_mineral_[posi] -= dt*sgn*fac1*fac2;
    }

    // Finally, update mineral concentrations on the "outside"
    for(int i = 0; i < ICS_->size_sup_min_; ++i)
    {
        const double mineral_conc = ctot_mineral_[ICS_->pos_sup_min_[i]];
        c_min[ICS_->krmap_[i]]    =  mineral_conc;
    }
}

void BasVec::update_mineral_concentrations_equilibrium(const double* c0_min,
                                                       const std::array<double, 3>& mass_phase,
                                                       double* c_min)
{
    for (int i = 0; i < ICS_->size_rock_; ++i)
    {
        const int mineral_index = ICS_->krmap_[ICS_->size_sup_min_ + i];
        const double dc_min = ICS_->SM_mineral_->delta_[ICS_->pos_buffer_[i]]
                              * mass_phase[gFluidPhase::WATER]
                              / mass_phase[ICS_->phase_[i]];

        c_min[mineral_index] = c0_min[mineral_index] + dc_min;
    }
}

void BasVec::update_ioxch_complex()
{
    if(!ICS_->includes_ion_exchange()) return;

    for (int i = 0; i < static_cast<int>(pos_exch_comp_.size()); ++i)
    {
        const double QK = ctot_[pos_exch_bas_[i]]/ fabs(ICS_->SM_all_->charge_[pos_exch_comp_[i]]);
        ICS_->SM_all_->log_QK_[pos_exch_comp_[i]] = std::log10(QK);
    }
}

void BasVec::update_dG(hkf& HKF_EOS, double T, double P)
{
    HKF_EOS.WaterProp(T, P);
    ICS_->CP_.ew_ = HKF_EOS.epsw_;
    ICS_->CP_.rho_w_ = HKF_EOS.rhow_;

    ChemTable& SM_mineral = *ICS_->SM_mineral_.get();
    const int size_hkf_min = SM_mineral.noRows_ - SM_mineral.size_analytical_;
    HKF_EOS.dGMineral(T,
                      P,
                      SM_mineral.deltaG_.data(),
                      SM_mineral.S_.data(),
                      SM_mineral.a1_.data(),
                      SM_mineral.a2_.data(),
                      SM_mineral.a3_.data(),
                      SM_mineral.mol_volume_.data(),
                      size_hkf_min,
                      SM_mineral.dG_TP_.data()
    );

    ChemTable& SM_basis = *ICS_->SM_basis_.get();
    const int size_hkf_basis = SM_basis.noRows_ - SM_basis.size_analytical_;
    HKF_EOS.dGIons(T,
                   P,
                   SM_basis.deltaG_.data(),
                   SM_basis.S_.data(),
                   SM_basis.a1_.data(),
                   SM_basis.a2_.data(),
                   SM_basis.a3_.data(),
                   SM_basis.a4_.data(),
                   SM_basis.c1_.data(),
                   SM_basis.c2_.data(),
                   SM_basis.omega_.data(),
                   SM_basis.charge_.data(),
                   SM_basis.re_.data(),
                   size_hkf_basis,
                   pos_water_,
                   SM_basis.dG_TP_.data(),
                   SM_basis.mol_volume_.data()
    );
    SM_basis.dG_TP_[pos_water_] = HKF_EOS.G_;  // H2O

    ChemTable& SM_all = *ICS_->SM_all_.get();
    const int size_hkf_all = SM_all.noRows_ - SM_all.size_analytical_;
    HKF_EOS.dGIons(T,
                   P,
                   SM_all.deltaG_.data(),
                   SM_all.S_.data(),
                   SM_all.a1_.data(),
                   SM_all.a2_.data(),
                   SM_all.a3_.data(),
                   SM_all.a4_.data(),
                   SM_all.c1_.data(),
                   SM_all.c2_.data(),
                   SM_all.omega_.data(),
                   SM_all.charge_.data(),
                   SM_all.re_.data(),
                   size_hkf_all,
                   -1,
                   SM_all.dG_TP_.data(),
                   SM_all.mol_volume_.data()
    );
}

void BasVec::update_logK([[maybe_unused]] hkf& HKF_EOS, double T)
{
    const double fact = 1.0 / (NumericalConstants::LNTEN*PhysicalConstants::IdealGasConstant*T);
    const auto& dG_basis = ICS_->SM_basis_->dG_TP_;

    ChemTable& SM_mineral = *ICS_->SM_mineral_.get();
    const int size_hkf_min = SM_mineral.noRows_ - SM_mineral.size_analytical_;
    const auto& dG_min = ICS_->SM_mineral_->dG_TP_;
    auto& logK_min = ICS_->SM_mineral_->logK_;

    for (int i = 0; i < size_hkf_min; ++i)
    {
        logK_min[i] = dG_min[i];
        for (int j = 0; j < SM_mineral.noColumns_; ++j)
        {
            logK_min[i] -= SM_mineral.M_[i][j] * dG_basis[j];
        }
        logK_min[i] *= fact;
    }

    const int pos_water_vapor = SM_mineral.get_row_index("H2O,g");
    if (pos_water_vapor >= 0)
    {
        const double water_vapor_logK = HKF_EOS.waterVaporLogK(T, PhysicalConstants::standard_gas_pressure);
        logK_min[pos_water_vapor] = water_vapor_logK;
        SM_mineral.dG_TP_[pos_water_vapor] = dG_basis[pos_water_] + water_vapor_logK / fact;
    }
    SM_mineral.update_logK_analytical(T);

    // Other species (ions, aqueous complexes, surface complexes, ...)
    ChemTable& SM_all = *ICS_->SM_all_.get();
    const int size_hkf_all = SM_all.noRows_ - SM_all.size_analytical_;
    const auto& dG_all = ICS_->SM_all_->dG_TP_;
    auto& logK_all = ICS_->SM_all_->logK_;

    for (int i = 0; i < size_hkf_all; ++i)
    {
        logK_all[i] = dG_all[i];
        for (int j = 0; j < SM_all.noColumns_; ++j)
        {
            logK_all[i] -= SM_all.M_[i][j] * dG_basis[j];
        }
        logK_all[i] *= fact;
    }
    SM_all.update_logK_analytical(T);
}
// Assume equilibrium phases behaves as ideal gas and calculate molar volume at reference pressure
void BasVec::update_molar_volume_gases_equilibrium_phases()
{
//    int pos = ICS_->SM_mineral_->
    int pos = ICS_->SM_mineral_->get_row_index("CO2,g");
    if(pos > -1)
        ICS_->SM_mineral_->mol_volume_[pos] = PhysicalConstants::IdealGasConstant * PhysicalConstants::ambient_temperature / PhysicalConstants::standard_gas_pressure;
}

void BasVec::set_temperature_for_equilibrium_reactions(hkf& HKF_EOS)
{
    const double ti = Temp_[gFluidPhase::WATER];
    const double pi = Pres_[gFluidPhase::WATER];

    ChemTable& SM_mineral = *ICS_->SM_mineral_;
    ChemTable& SM_all = *ICS_->SM_all_;

    SM_mineral.T_curr_ = ti;
    SM_mineral.P_curr_ = pi;
    SM_all.T_curr_ = ti;
    SM_all.P_curr_ = pi;



    update_dG(HKF_EOS, ti, pi);
    update_logK(HKF_EOS, ti);

    // Debye-Huckel
    static constexpr double CAdh[] = { 0.51, -1.154, 35.697, -182.023, 346.528 };
    static constexpr double CBdh[] = { 0.325, 0.08, 1.441, -6.541, 11.655 };

    ICS_->CP_.Adh_= 0.0;
    ICS_->CP_.Bdh_ = 0.0;
    for(int j=0; j < 5; ++j)
    {
        ICS_->CP_.Adh_ += pow(0.001*(ti-273.15), j)*CAdh[j];
        ICS_->CP_.Bdh_ += pow(0.001*(ti-273.15), j)*CBdh[j];
    }
    ICS_->CP_.F_div_RT_ = PhysicalConstants::Faraday/(PhysicalConstants::IdealGasConstant*ti);
}
// used in the initialization only, and the ideal gas law is assumed
void BasVec::add_gas_phase_to_total_conc()
{
    if (size_gase_phase_ < 1) return;
    [[maybe_unused]] static constexpr double Rg = PhysicalConstants::IdealGasConstant;
    const double fact = 1e-3*gas_phase_frac_ / (PhysicalConstants::IdealGasConstant * Temp_[gFluidPhase::GAS]); //mol/L
    for (int i = 0; i < size_; ++i)
    {
        for (int j = 0; j < size_gase_phase_; ++j)
        {
            if (i != pos_exp_ && i != pos_pH_ && i != pos_pe_ && i != pos_water_)
            {
                ctot_[i] += ICS_->SM_mineral_->M_[pos_gas_phase_[j]][i] * gas_phase_pressure_[j] * fact;
            }
        }

    }
}
void BasVec::set_temperature_for_mineral_kinetics()
{
    if(size_sup_min_ == 0) return;

    const double dT = 1.0 / Temp_[gFluidPhase::WATER] - 1.0 / 298.15;
    static constexpr double Rg = PhysicalConstants::IdealGasConstant;

    for(int i=0; i < size_sup_min_; ++i)
    {
        const double k1_Tref = ICS_->rate_[pos_rel_sup_min_[i]]->k1_;
        const double k2_Tref = ICS_->rate_[pos_rel_sup_min_[i]]->k2_;
        const double Ea1 = ICS_->rate_[pos_rel_sup_min_[i]]->Ea1_;
        const double Ea2 = ICS_->rate_[pos_rel_sup_min_[i]]->Ea2_;
        const double SA = ICS_->rate_[pos_rel_sup_min_[i]]->SA_;
        rate_laws_[i].k1_ = SA*k1_Tref * exp(-Ea1 / Rg * dT);
        rate_laws_[i].k2_ = SA*k2_Tref * exp(-Ea2 / Rg * dT);
    }
}

/** Calculates activities based on the Debye-Hückel formulation. */
void BasVec::update_activity_coefficients(double Io)
{
    const double Io_sqrt = sqrt(Io);
    const double A = ICS_->CP_.Adh_;
    const double B = ICS_->CP_.Bdh_;
    const double  Setschenow_b = 0.1;

    const ChemTable& SM_basis = *ICS_->SM_basis_;
    for (int i = 0; i < size_; ++i)
    {
        if(i != pos_pe_)
            log_g_[i] = -A * SM_basis.charge_[i] * SM_basis.charge_[i] * Io_sqrt
                        / (1.0 + B * SM_basis.a0_[i] * Io_sqrt);
        if (SM_basis.charge_[i] == 0.)
        {
            log_g_[i] = Setschenow_b * Io;
        }
    }

    ChemTable& SM_all = *ICS_->SM_all_;
    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        SM_all.log_g_[i] = -A * SM_all.charge_[i] * SM_all.charge_[i] * Io_sqrt
                           / (1.0 + B * SM_all.a0_[i] * Io_sqrt);
        if (SM_all.charge_[i] == 0.)
        {
            SM_all.log_g_[i] = Setschenow_b * Io;
        }
    }
}

bool BasVec::new_temperature() const
{
    return ICS_->SM_all_->T_curr_ != Temp_[gFluidPhase::WATER];
}

bool BasVec::new_pressure() const
{
    return ICS_->SM_all_->P_curr_ != Pres_[gFluidPhase::WATER];
}

double BasVec::calc_pH() const
{
    assert(pos_pH_ > -1);
    return -log_a_[pos_pH_];
}

double BasVec::calc_exchanger_charge() const
{
    // Note: basis species NOT included (e.g., X-)
    const ChemTable& SM_all = *ICS_->SM_all_;
    double cbal_io = 0.0;
    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if(SM_all.type_[i] == GeochemicalComponentType::ION_EXCHANGE)
        {
            cbal_io += SM_all.charge_[i]*POW10(SM_all.log_m_[i]);
        }
    }
    return cbal_io;
}

double BasVec::calc_solution_charge() const{

    const ChemTable& SM_basis = *ICS_->SM_basis_;
    const ChemTable& SM_all = *ICS_->SM_all_;

    double cbal_aq = 0.0;
    for (int i = 0; i < size_; ++i)
    {
        if(SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            cbal_aq += SM_basis.charge_[i]*POW10(log_m_[i]);
        }
    }

    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if(SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            cbal_aq += SM_all.charge_[i]*POW10(SM_all.log_m_[i]);
        }
    }
    return cbal_aq;
}

// calculates partial pressures of a mixture of gasses,
// currently:: assumes gases behaves ideally and ignore capillary pressure
void BasVec::calculate_gas_phase_pressure()
{
    double Pres_gas = 0.;
    for (int i = 0; i < size_gase_phase_; ++i)
    {
        int pos = pos_gas_phase_[i];
        Pres_gas += POW10(ICS_->SM_mineral_->log_a_[pos]) / ICS_->SM_mineral_->fugacity_[pos];
    }
    Pres_gas *= PhysicalConstants::standard_gas_pressure;
    Pres_[gFluidPhase::GAS] = Pres_[gFluidPhase::OIL] = Pres_[gFluidPhase::WATER] = Pres_gas;
}

double BasVec::calc_surface_charge(bool coulombs_per_square_meter) const{

    const ChemTable& SM_all = *ICS_->SM_all_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    double surface_charge = 0.0;
    for (int i = 0; i < size_; ++i)
    {
        surface_charge += SM_basis.scharge_[i]*POW10(log_m_[i]);  // basis species
    }
    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        surface_charge += SM_all.scharge_[i]*POW10(SM_all.log_m_[i]);  // surface complexes
    }

    if(coulombs_per_square_meter)
    {
        return surface_charge*(PhysicalConstants::Faraday / SA_);
    }
    return surface_charge;
}

std::size_t BasVec::indexOfBasisSurfaceComplex(std::size_t indexOfSecondarySurfaceComplex) const
{
    const ChemTable& SM_all = *ICS_->SM_all_;
    const ChemTable& SM_basis = *ICS_->SM_basis_;

    const auto& complex_name = SM_all.row_name_[indexOfSecondarySurfaceComplex];

    bool found_basis_complex = false;
    for (int i = 0; i < SM_basis.noRows_; ++i)
    {
        const auto stoichiometric_coeff = SM_all.M_[indexOfSecondarySurfaceComplex][i];
        if (SM_basis.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX && stoichiometric_coeff != 0 && SM_basis.row_name_[i] != "E")
        {
            if (found_basis_complex)
            {
                const std::string error_msg = "The secondary surface complex " + complex_name + " can only be associated with a single basis surface complex...";
                throw InvalidInputException(error_msg);
            }
            else return i;
        }
    }

    if (!found_basis_complex)
    {
        const std::string error_msg = "Error: The secondary surface complex " + complex_name + " has no associated basis surface complex...";
        throw InvalidInputException(error_msg);
    }
    return -1;
}

/**
 * This routine assign one basis specie to each mineral phase, which presupposes
 * that all minerals are linearly independent. The trick is to find the mineral
 * with the least number of species first, and run through the list in that order.
*/
void BasVec::assign_basis2mineral()
{
    static constexpr int MAX_SPECIES_IN_MINERAL = 1000;

    std::vector<int> ignore_list;
    std::vector<int> ignore_list_buf;
    ignore_list.reserve(size_);
    ignore_list_buf.reserve(size_);

    // The following basis species will never be coupled to minerals for the purpose of basis switching.
    //
    if(charge_balance_)
        ignore_list.push_back(pos_pH_);  // H+ (Before 31/5-2022: Only if charge_balance_)
    ignore_list.push_back(pos_water_);  // H2O
    if (pos_exp_ > -1) ignore_list.push_back(pos_exp_);  // Special species E.
 //   if (pos_pe_ > -1) ignore_list.push_back(pos_pe_);  // Special species e- (electron).

    for(int mineral_count = 0; mineral_count < ICS_->size_min_; ++mineral_count)
    {
        int no_basis_species_left = MAX_SPECIES_IN_MINERAL;
        int mineral_index = 0;
        int pos_i = -1;
        int pos_b = -1;
        int pos_bas = -1;

        for (int i = 0; i < ICS_->size_min_; ++i)
        {
            int no_basis_species_in_mineral;  // (Other than those we ignore)
            if(element_is_in_container(ICS_->pos_min_[i], ignore_list_buf))
            {
                no_basis_species_in_mineral = MAX_SPECIES_IN_MINERAL;
            }
            else
            {
                no_basis_species_in_mineral = ICS_->SM_mineral_->find_no_basis_species(ICS_->pos_min_[i],
                                                                                       &pos_bas,
                                                                                       ignore_list);
            }

            if (no_basis_species_in_mineral > 0 && no_basis_species_in_mineral < no_basis_species_left)
            {
                pos_i                   = pos_bas;
                mineral_index           = i;
                no_basis_species_left   = no_basis_species_in_mineral;
                pos_b                   = ICS_->pos_min_[i];
            }
        }
        assert(pos_i > -1);
        assert(pos_b > -1);

        if (ICS_->PRINT_DEBUG_CHEM_ >= DebugInfoLevel::ALOT)
        {
            fmt::print(
                "Number of basis species left in mineral {} is {}, choose basis species {}.",
                ICS_->SM_mineral_->row_name_[ICS_->pos_min_[mineral_index]],
                no_basis_species_left,
                ICS_->SM_mineral_->col_name_[pos_i]
            );
        }

        ICS_->pos_min_bas_[mineral_index] = pos_i;
        ignore_list.push_back(pos_i);  // Ignore the basis specie from here on
        ignore_list_buf.push_back(pos_b);  // Also ignore buffer mineral
    }

    for (int i = 0; i < ICS_->size_sup_min_; ++i)
    {
        ICS_->pos_sup_bas_[i] = ICS_->pos_min_bas_[i];
    }
    for (int i = 0; i < ICS_->size_rock_;++i)
    {
        ICS_->pos_rock_[i] = ICS_->pos_min_bas_[ICS_->size_sup_min_+i];
    }
}

// Note that sm_min_ is created separately...
void BasVec::update_basis_buffer_matrices(int old_size_rock)
{
    // Note that we always create this matrix, even with no basis-switching.
    // Do we need do to that?
    if (!beta_bas_inv_)
    {
        // Create basis transformation matrix
        beta_bas_inv_ = allocateMemoryForMatrix(size_);
        fillIdentityMatrix(beta_bas_inv_, size_);
    }

    assert(size_rock_ >= 0);
    if(size_rock_ == 0)
    {
        assert(old_size_rock == 0);  // We never remove buffers...
        return;
    }

    if(old_size_rock > 0)
    {
        assert(old_size_rock == size_rock_ - 1);  // We never add more than a single new buffer at a time.
        freeMatrixMemory(sm_buf_, old_size_rock);  // Old matrix has wrong dimensions.
    }

    update_beta_inv();
    sm_buf_ = allocateMemoryForMatrix(size_rock_);
    update_sm_buf(size_rock_, pos_buffer_.data(), pos_rock_.data(), sm_buf_);
}

void BasVec::update_sm_buf(int size_r, int* pos_b, int* pos_r, double** sm)
{
    ChemTable& SM_mineral = *ICS_->SM_mineral_;
    SM_mineral.in_fmatrix(pos_b, size_r, pos_r, size_r, sm);

    // We need the inverse of the transposed reduced stoichiometric matrix
    double** m_dum = allocateMemoryForMatrix(size_r);
    for(int i=0; i < size_r; ++i)
    {
        for(int j=0; j < size_r; ++j)
        {
            m_dum[i][j] = sm[j][i];
        }
    }

    double determinant = -1;
    invert_matrix(m_dum, size_r, sm, determinant);
    freeMatrixMemory(m_dum, size_r);
}

/*
* Finds the basis transformation matrix, calculates its inverse, and stores it in beta_inv.
* Prior to calling this function, it is assumed that memory space for beta_bas_inv has been allocated.
*/
void BasVec::update_beta_inv()
{
    assert(beta_bas_inv_ != nullptr);  // Throw an error if this condition somehow should fail?

    const ChemTable& SM_mineral = *ICS_->SM_mineral_;

    // Allocate temporary workspace (deleted inside this function)
    double** beta_bas = allocateMemoryForMatrix(size_);

    // Basis transformation matrix
    for(int j=0; j < size_; ++j)
    {
        // Check if basis species j is one of the rock buffered species.
        int length = index_of_element_in_array(j, pos_rock_.data(), size_rock_);

        if(length >= 0)
        {
            for(int i=0; i < size_; ++i)
            {
                beta_bas[j][i] = SM_mineral.M_[pos_buffer_[length]][i];
            }
        }
        else
        {
            for(int i=0; i < size_; ++i)
            {
                if(i==j)    beta_bas[j][i] = 1.0;
                else        beta_bas[j][i] = 0.0;
            }
        }
    }


    double determinant;
    invert_matrix(beta_bas, size_, beta_bas_inv_, determinant);

    // Calculate transformed stoichiometric matrix for complexes
    if (M_all_trans_.empty())
    {
        M_all_trans_.resize(ICS_->SM_all_->noRows_);
        for (int i = 0; i < ICS_->SM_all_->noRows_; ++i)
        {
            M_all_trans_[i].resize(size_);
        }
    }

    for (int i = 0; i < ICS_->SM_all_->noRows_; ++i)
    {
        for (int j = 0; j < size_; ++j)
        {
            M_all_trans_[i][j] = 0.0;

            // Basis transformation
            for (int p = 0; p < size_; ++p)
            {
                M_all_trans_[i][j] += ICS_->SM_all_->M_[i][p] * beta_bas_inv_[p][j];
            }
        }
    }

    freeMatrixMemory(beta_bas, size_);
}

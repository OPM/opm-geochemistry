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
#include <opm/simulators/geochemistry/Core/ChemGCSolver.h>

GCSolver::GCSolver(GeoLogger* logger, BasVec* Vchem, int sizeB, int sizeM, int MaxSizeBasis, int MaxSizeComp, int MaxSizeMin)
    : c0_()
    , c0_ads_()
    , c0_min_()
    , cT_()
    , cT_ads_()
    , cT_ads_io_()
    , cT_ads_sc_()
    , cT_DL_()
    , cT_min_()
    , SI_min_()
    , mass_phase_{ 1.0, 1.0, 1.0 }
    , sizeB_(sizeB)
    , sizeM_(sizeM)
    , MaxSizeBasis_(MaxSizeBasis)
    , MaxSizeComp_(MaxSizeComp)
    , MaxSizeMin_(MaxSizeMin)
    , options_()
    , Vchem_(Vchem)
    , HKF_EOS_(nullptr)
    , logger_(logger)
//
    , c_dl_excess_()
{
    if(!logger_)
    {
        error_and_exit("ChemGCSolver must be provided with a Logger...");
    }

    HKF_EOS_ = std::make_unique<hkf>();

    c0_.resize(sizeB_);
    c0_ads_.resize(sizeB_);
    c0_min_.resize(sizeM_);

    cT_.resize(sizeB_);
    cT_DL_.resize(sizeB_);
    cT_ads_.resize(sizeB_);
    cT_ads_io_.resize(sizeB_);
    cT_ads_sc_.resize(sizeB_);
    cT_min_.resize(sizeM_);
    SI_min_.resize(sizeM_);

    alpha_a_gas_phase_.resize(sizeB_);


    smpow_.resize(MaxSizeComp + 1);
    mpow_.resize(MaxSizeBasis + 1);
    Fchem_.resize(MaxSizeBasis + 1);
    x0_.resize(MaxSizeBasis + 1);
    x_init_.resize(MaxSizeBasis + 1);
    delta_calc_.resize(MaxSizeBasis + 1);
    chem_indx_.resize(MaxSizeBasis + 1);
    c_dl_excess_.resize(MaxSizeBasis + 1);

    Jacobi_ = allocateMemoryForMatrix(MaxSizeBasis+1);
    Jacobi_dum_ = allocateMemoryForMatrix(MaxSizeBasis+1);

    wksp_linalg_.resize(MaxSizeBasis+1);
}

GCSolver::~GCSolver()
{
    if(Jacobi_) freeMatrixMemory(Jacobi_, MaxSizeBasis_ + 1);
    if(Jacobi_dum_) freeMatrixMemory(Jacobi_dum_, MaxSizeBasis_ + 1);
}

bool GCSolver::massBalanceHasConverged() const {
    return solver_state_.massBalanceStatus_ == ConvergenceStatus::CONVERGED;
}

void GCSolver::write_convergence_status(const std::string& name,bool converged )
{
    std::ofstream file(name);
    if (file.is_open())
    {
        if (converged)
            file << 1 << std::endl;
        else
            file << 0 << std::endl;
        file.close();
    } else {
     std::cerr << "Error writing convergence file!" << std::endl;
    }
}

/** Note: ordering is the same as in ICS_full_. */
bool GCSolver::SolveChem()
{
    const InitChem& ICSp = *Vchem_->ICS_;

    for (int i = 0; i < Vchem_->size_gas_; ++i)
    {   // Currently, no fugacity model has been implemented.
        ICSp.SM_mineral_->fugacity_[Vchem_->pos_gas_[i]] = 1.0;
    }
    // Gas phase
    // maybe better to do in the initialization of Vchem?
    // @ah Vchem_->add_gas_phase_to_total_conc();

 //   Vchem_->update_molar_volume_gases_equilibrium_phases();

    // Update logK-values and rate constants

    Vchem_->set_temperature_for_equilibrium_reactions(*HKF_EOS_);
    Vchem_->set_temperature_for_mineral_kinetics();




    Vchem_->update_ioxch_complex();

    bool converged=solve_chemistry();

    for (int i = 0; i < ICSp.size_basis_; ++i)
    {
        const int speciesIdx = ICSp.kcmap_[i];
        const int basVecIdx = ICSp.pos_[i];

        cT_[speciesIdx] = Vchem_->ctot_calc_[0][basVecIdx];
        cT_ads_sc_[speciesIdx] = Vchem_->ctot_calc_[1][basVecIdx];
        cT_ads_io_[speciesIdx] = Vchem_->ctot_calc_[2][basVecIdx];
        cT_ads_[speciesIdx] = cT_ads_sc_[speciesIdx] + cT_ads_io_[speciesIdx];
    }
    return converged;
}

/**
 * IMPORTANT NOTE (1/12-22):
 *    The signature of this function has changed, c_bas is now constant
 *    and only refers to the basis concentrations.
 *
 * It is the responsibility of the caller to ensure that all of the
 * pointers point an arrays / vector with the right size.
 *
 * @param [in] key Old mineral key, indicating which minerals are present.
 * @param [in] dt Time step used.
 * @param [in] mass_phase Mass of, respectively, the water, oil, and gas phases.
 * @param [in] c_bas Old total concentrations of basis species.
 * @param [in] c_ads Old adsorbed concentrations (only used if relevant).
 * @param [in,out] c_min Mineral concentrations (old --> new).
 * @param [in,out] log_a_min  Log activities for minerals (old --> new).
 * @param [in,out] output_values At exit, contains concentration differences (old-new),
 *                               as well as surface charge and surface potential (when relevant),
 *                               plus pH.
 *
 * @return Potentially new mineral key (or the old one back).
 */
int GCSolver::SolveChem_ST(int key,
                           [[maybe_unused]] double dt,
                           const std::array<double, 3>& mass_phase,
                           const double* c_bas,
                           const double* c_ads,
                           double* c_min,
                           double* log_a_min,
                           double* output_values)
{
    const ChemTable& SM_mineral = *Vchem_->ICS_->SM_mineral_;
    InitChem* const ICSp = Vchem_->ICS_;

    mass_phase_ = mass_phase;

    for (int i = 0; i < sizeB_; ++i)
    {
        cT_[i] = c_bas[i];
        c0_[i] = c_bas[i];

        output_values[i] = 0.0;
        if (ICSp->includes_ion_exchange() || ICSp->includes_surface_complexes())
        {
            cT_ads_[i]          = c_ads[i];
            c0_ads_[i]          = c_ads[i];
            output_values[sizeB_ + i]   = 0.0;
        }
    }

    for (int i = 0; i < sizeM_; ++i)
    {
        cT_min_[i] = c_min[i];
        c0_min_[i] = c_min[i];
    }

    // c_bas and c_min are sorted as in ICS_full
    for (int i = 0; i < ICSp->size_basis_; ++i)
    {
        const int posi = ICSp->pos_[i];
        const int posi_kcmap = ICSp->kcmap_[i];

        const double ct = std::max(NumericalConstants::ZERO_THRESHOLD, c_bas[posi_kcmap]);
        Vchem_->ctot_[posi] = ct;
        ICSp->c_vchem_[i] = ct;

        if (posi != Vchem_->pos_pH_)
        {
            // For species other than H+, use the total conc. as initial guess for the basis conc.
            const double log_mi = log10(Vchem_->ctot_[posi]);
            Vchem_->log_m_[posi] = log_mi;
            Vchem_->log_a_[posi] = log_mi;
        }
    }

    Vchem_->update_ioxch_complex();

    for (int i = 0; i < ICSp->size_sup_min_; ++i)
    {
        Vchem_->ctot_mineral_[ICSp->pos_sup_min_[i]] = ICSp->c_sup_min_[i] = c_min[ICSp->krmap_[i]];
        ICSp->SM_mineral_->log_af_[ICSp->pos_sup_min_[i]] = log_a_min[ICSp->krmap_[i]];
    }

    for (int i = 0; i < ICSp->size_rock_; ++i)
    {
        const double buffer_conc = c_min[ICSp->krmap_[i + ICSp->size_sup_min_]];
        const double buffer_log_a_min = log_a_min[ICSp->krmap_[i + ICSp->size_sup_min_]];

        Vchem_->ctot_mineral_[ICSp->pos_buffer_[i]]         = buffer_conc;
        ICSp->c_buffer_[i]                                  = buffer_conc;
        ICSp->SM_mineral_->log_af_[ICSp->pos_buffer_[i]]    = buffer_log_a_min;
        ICSp->SM_mineral_->delta_[ICSp->pos_buffer_[i]]     = 0.0;
    }

    for (int i = 0; i < Vchem_->size_gas_; ++i)
    {
        ICSp->SM_mineral_->fugacity_[Vchem_->pos_gas_[i]]   = 1.0;  // Not implemented.
    }

    if(Vchem_->new_temperature() || Vchem_->new_pressure())
    {
       Vchem_->set_temperature_for_equilibrium_reactions(*HKF_EOS_);
    }

    // The number of rate-minerals can change even at constant temperature.
    Vchem_->set_temperature_for_mineral_kinetics();

    const bool chemistry_converged = solve_chemistry();
    if (!chemistry_converged)
    {
        return key;  // return old key in case of failure
    }

    if (ICSp->PRINT_DEBUG_CHEM_ == DebugInfoLevel::MAX)
    {
        // We print out a crazy number of files...
        writeCurrentBasVecToFile("sol_I" + std::to_string(solver_state_.noCalls_) + ".out");
    }

    // Non-linear rate equations
    // 09.05.2022: New key should be generated for any mineral type
    //             --> Changed all instances of size_sup_min_ to size_min_
    //
    if (!Vchem_->equilibrate_ && ICSp->size_min_ > 0 && (Vchem_->size_min_ != ICSp->size_min_))
    {
        const int mineralIndex = Vchem_->index_of_most_supersaturated_mineral(ICSp->pos_min_,
                                                                              ICSp->size_min_,
                                                                              Vchem_->pos_min_,
                                                                              Vchem_->size_min_);
        if (mineralIndex >= 0)
        {
            logger_->info("Adding buffer mineral {}...", SM_mineral.row_name_[ICSp->pos_min_[mineralIndex]]);
            if (ICSp->PRINT_DEBUG_CHEM_) Vchem_->write("debugi.out");

            const int mineralIndex_relative_to_ICS_full = ICSp->krmap_[mineralIndex];
            c_min[mineralIndex_relative_to_ICS_full] = NumericalConstants::ZERO_THRESHOLD;

            return create_unique_key_from_positive_elements(c_min, sizeM_);
        }
    }

    // If no new key was generated, we proceed...

    // Calculate new mineral concentrations...
    Vchem_->update_mineral_concentrations_kinetics(dt_, c_min);
    Vchem_->update_mineral_concentrations_equilibrium(c0_min_.data(), mass_phase_, c_min);
    for(int i=0; i < sizeM_; ++i)
    {
        cT_min_[i] = c_min[i];
    }

    for (int i = 0; i < Vchem_->size_gas_; ++i)
    {
        log_a_min[ICSp->krmap_[Vchem_->pos_rel_gas_[i] + ICSp->size_sup_min_]] = ICSp->SM_mineral_->log_af_[Vchem_->pos_gas_[i]];
    }

    for (int i = 0; i < ICSp->size_sup_min_; ++i)
    {
        SI_min_[ICSp->krmap_[i]] = ICSp->SM_mineral_->log_a_[ICSp->pos_sup_min_[i]];
    }

    for (int i = 0; i < ICSp->size_rock_; ++i)
    {
        SI_min_[ICSp->krmap_[ICSp->size_sup_min_ + i]] = ICSp->SM_mineral_->log_a_[ICSp->pos_buffer_[i]];
    }

    for (int i = 0; i < ICSp->size_basis_; ++i)
    {
        const double frac_DL = diffuse_layer_props_.frac_DL_;

        const double c_free = Vchem_->ctot_calc_[0][ICSp->pos_[i]];
        const double c_sc = Vchem_->ctot_calc_[1][ICSp->pos_[i]];
        const double c_io = Vchem_->ctot_calc_[2][ICSp->pos_[i]];
        const double cdl_excess = Vchem_->ctot_dl_excess_[ICSp->pos_[i]];

        const int speciesIdx = ICSp->kcmap_[i];
        cT_DL_[speciesIdx] = cdl_excess + frac_DL * c_free;  // Before: c_free + cdl_excess
        cT_ads_sc_[speciesIdx] = c_sc;
        cT_ads_io_[speciesIdx] = c_io;
        cT_ads_[speciesIdx] = cT_ads_sc_[speciesIdx] + cT_ads_io_[speciesIdx] + cT_DL_[speciesIdx];  // Before: c_sc + c_io + f_DL*cT_DL
        cT_[speciesIdx] = (1.0 - frac_DL)*c_free + cT_ads_[speciesIdx];

        // NB: At exit, we store concentration DIFFERENCES
        output_values[speciesIdx] = c0_[speciesIdx] - cT_[speciesIdx];

        if (ICSp->includes_ion_exchange() || ICSp->includes_surface_complexes())
        {
            // Note: The first time we enter the for loop, we overwrite dt (or
            //       0.0, if no dt was input). Thus, the indices of adsorbed
            //       species are all shifted to the left by 1, compared to when
            //       the function was called...
            output_values[sizeB_ + speciesIdx] = c_ads[speciesIdx] - cT_ads_[speciesIdx];
        }
    }

    const double pH = -Vchem_->log_a_[Vchem_->pos_pH_];
    if (ICSp->includes_surface_complexes())
    {
        output_values[2 * sizeB_]       = Vchem_->sch_;
        output_values[2 * sizeB_ + 1]   = Vchem_->psi_;
        output_values[2 * sizeB_ + 2]   = pH;
    }
    else if (ICSp->includes_ion_exchange())     output_values[2 * sizeB_]   = pH;
    else                                        output_values[sizeB_]       = pH;

    return key; // Return old key
}

bool GCSolver::solve_chemistry()
{
    ++solver_state_.noCalls_;

    ChemTable& SM_mineral = *Vchem_->ICS_->SM_mineral_;
    ChemTable& SM_all = *Vchem_->ICS_->SM_all_;

    solver_state_.reset_pH();

    if (Vchem_->surface_options_.CALCULATE_DIFFUSE_LAYER)
    {
        fill_zero(Vchem_->ctot_dl_excess_);
        fill_zero(c_dl_excess_);
    }

    // Initialize aqueous complexes.
    // This is necessary in cases where no minerals are present to buffer the pH.
    // Small changes in initial conditions will make pH oscillate slightly.
    // This is expected from a physical point of view...
    fill_value(SM_all.log_a_, NumericalConstants::LOG10_ALMOST_ZERO);
    fill_value(SM_all.log_m_, NumericalConstants::LOG10_ALMOST_ZERO);
    fill_zero(SM_all.log_g_);
    fill_zero(Vchem_->log_g_);


    // =========================== MAIN SOLVER (pH) LOOP ========================
    using SecantMethod = SecantMethodReal<double, DoNothingOnError>;

    double Io_old = Vchem_->update_ionic_strength();
    Vchem_->update_activity_coefficients(Io_old);
    double Z_factor_old = Z_factor_gas_;
    if (Vchem_->surface_options_.CALCULATE_DIFFUSE_LAYER) integrate_diffuse_layer_for_mono_or_divalent_ions();

    // Experimental model with variable cation / anion - exchange capacity
    if (Vchem_->surface_options_.EXCHANGE_MODEL_DL)
    {
        Vchem_->calc_ctot_aq(); // needed to estimate exchange capacities
        Vchem_->update_exchange_capacity_dl_model(diffuse_layer_props_.frac_DL_);
        Vchem_->update_ioxch_complex(); // update logQ_K on exchange species

        const int posXDL = Vchem_->ICS_->SM_basis_->get_row_index("XDL-");
        const int posYDL = Vchem_->ICS_->SM_basis_->get_row_index("YDL+");
        Vchem_->log_a_[posXDL] = log10(Vchem_->ctot_[posXDL]);
        Vchem_->log_m_[posXDL] = log10(Vchem_->ctot_[posXDL]);
        Vchem_->log_a_[posYDL] = log10(Vchem_->ctot_[posYDL]);
        Vchem_->log_m_[posYDL] = log10(Vchem_->ctot_[posYDL]);
    }

    double x_curr;  // == -pH
    double x_next{};
    double f_next;
    double x_old{};
    double f_old{};
    double delta_x;

    double exch_capacity_dl_old = 0.0;
    double exch_capacity_dl_new = 0.0;

    if (Vchem_->charge_balance_)
    {
        x_old = Vchem_->log_a_[Vchem_->pos_pH_];
        x_next = x_old * options_.PH_SECANT_METHOD_X1_FAC;
        f_old = f_pH_charge_balance(x_old);
        f_next = f_pH_charge_balance(x_next);
        delta_x = SecantMethod::computeStepSize(x_old, x_next, f_old, f_next);
        delta_x = restrictByMaximumAbsoluteValue(delta_x, options_.PH_MAX_STEP_);
    }
    else  // pH is fixed
    {
        const double negpH = Vchem_->log_a_[Vchem_->pos_pH_]; // -pH
        f_next = f_pH_charge_balance(negpH);
        x_curr = negpH;
    }

    double delta_dl = 0.0;
    // Helper lambda function used several times during pH iterations
    auto update_delta_dl_and_excess = [&]()
    {
      if (Vchem_->surface_options_.CALCULATE_DIFFUSE_LAYER)
      {
          calc_diffuse_layer_conc(c_dl_excess_.data());
          delta_dl = 0.0;
          for (int i = 0; i < Vchem_->size_; ++i)
          {
              delta_dl += std::fabs(Vchem_->ctot_dl_excess_[i] - c_dl_excess_[i]);
              Vchem_->ctot_dl_excess_[i] = c_dl_excess_[i];
          }
      }
    };

    bool has_converged = false;
    while(!has_converged)
    {
        ++solver_state_.noIterChargeBalance_;

        if (Vchem_->surface_options_.CALCULATE_DIFFUSE_LAYER)
        {
            integrate_diffuse_layer_for_mono_or_divalent_ions();
        }

        Io_old = Vchem_->Io_;
        const double Io_new = Vchem_->update_ionic_strength();
        Vchem_->update_activity_coefficients(Io_new);


        if (Vchem_->surface_options_.EXCHANGE_MODEL_DL)
        {
            Vchem_->calc_ctot_aq(); // needed to estimate exch.capacities
            exch_capacity_dl_old = exch_capacity_dl_new;
            exch_capacity_dl_new = Vchem_->update_exchange_capacity_dl_model(diffuse_layer_props_.frac_DL_);
            Vchem_->update_ioxch_complex(); // update logQ_K on exchange species
        }

        if (Vchem_->charge_balance_)
        {
            delta_x = SecantMethod::computeStepSize(x_old, x_next, f_old, f_next);
            delta_x = restrictByMaximumAbsoluteValue(delta_x, options_.PH_MAX_STEP_);
            const double x_tmp = x_next - delta_x;
            x_old = x_next;
            f_old = f_next;

            // Update solver state & check if ok
            x_next = x_tmp;
            f_next = f_pH_charge_balance(x_next);
            if (solver_state_.massBalanceStatus_ == ConvergenceStatus::NOT_CONVERGED)
            {
                return false;
            }
            update_delta_dl_and_excess();

            x_curr = x_old;
        }
        else
        {
            // Update solver state (incl. pH) & check if ok
            f_next = f_pH_charge_balance(Vchem_->log_a_[Vchem_->pos_pH_]);
            if (solver_state_.massBalanceStatus_ == ConvergenceStatus::NOT_CONVERGED)
            {
                return false;
            }
            update_delta_dl_and_excess();

            x_curr = Vchem_->log_a_[Vchem_->pos_pH_];
            f_next = 0.0; // since charge balance is not a constraint
        }

        //debug gas....

        if (Vchem_->size_gase_phase_ > 0 && fabs(f_next) < 1e-3)
        {
            Z_factor_old = Z_factor_gas_;
            //            update_fugacity_and_mol_volume();
           update_fugacity_and_mol_volume_mixtures();

           double alpha = 0.7;
           Z_factor_gas_ = alpha*Z_factor_old + (1-alpha)*Z_factor_gas_;
               // Vchem_->calculate_gas_phase_pressure();
           // if (Vchem_->new_pressure()) // recalculate pressure dependent logK values
           //{
           if(Vchem_->Pres_[gFluidPhase::WATER] < 100e6)
              Vchem_->set_temperature_for_equilibrium_reactions(*HKF_EOS_);
           //}

        }

        // Is this solution better than any solution we have seen before?
        const auto abs_fx = std::fabs(f_next);
        if (abs_fx < solver_state_.cached_best_f_pH_)
        {
            solver_state_.cached_best_pH_ = -x_curr;
            solver_state_.cached_best_f_pH_ = abs_fx;
        }

        const double delta_pH = Vchem_->charge_balance_ ? (x_curr - Vchem_->log_a_[Vchem_->pos_pH_]) : 0.0;
        const double delta_Io = Vchem_->Io_ - Io_old;
        const double PH_CONV_CRITERION = options_.PH_CONV_CRITERION_;
        const double DELTA_PH_CONV_CRITERION = options_.DELTA_PH_CONV_CRITERION;
        //const double charge_balance = Vchem_->charge_balance_ ? Vchem_->calc_solution_charge() : 0.0;
        double delta_exch_capacity_dl = exch_capacity_dl_new - exch_capacity_dl_old;

        has_converged =
            (fabs(delta_pH) < DELTA_PH_CONV_CRITERION) &&
            (fabs(delta_exch_capacity_dl) < PH_CONV_CRITERION) &&
            (abs_fx < PH_CONV_CRITERION) &&
            (fabs(delta_dl) < PH_CONV_CRITERION) &&
            //(fabs(charge_balance) < PH_CONV_CRITERION) &&
            (fabs(delta_Io) < PH_CONV_CRITERION &&
                fabs(Z_factor_old-Z_factor_gas_)<PH_CONV_CRITERION);

        logger_->info("Sch_={:12.8e}, psi_={:12.8e}, pH={:12.9e}, chrg={:12.9e}, Io_={}, Z= {:12.9e}, P(bar)= {:12.9e}",
                      Vchem_->sch_,
                      Vchem_->psi_,
                      Vchem_->calc_pH(),
                      f_next,
                      Vchem_->Io_,
                      Z_factor_gas_,
                       Vchem_->Pres_[gFluidPhase::GAS]*1e-5);

        if(has_converged)
        {
            solver_state_.massBalanceStatus_ = ConvergenceStatus::CONVERGED;
            break;
        }
        else if (solver_state_.noIterChargeBalance_ > options_.CHEM_MAX_PH_ITER_)
        {
            break;
        }
    }

    if (!has_converged)
    {
        std::string warning_msg("Solver did not converge within the set tolerance...");
        if (Vchem_->size_min_ == 0)
        {
            warning_msg += " Note that no minerals are present to buffer the pH.";
        }
        logger_->warning(warning_msg.c_str());

        if (Vchem_->ICS_->PRINT_DEBUG_CHEM_ == DebugInfoLevel::ALOT)
        {
            writeCurrentBasVecToFile("pH_dd" + std::to_string(solver_state_.noCalls_) + ".out");
        }

        logger_->info("Rerunning Newton-iterations a final time using the \"best\" pH...");
        f_next = f_pH_charge_balance(-solver_state_.cached_best_pH_);
        update_delta_dl_and_excess();
    }

    // POST-PROCESSING & OUTPUT
    Vchem_->calc_complex_sparse(SM_mineral);  // Saturation index of minerals
    Vchem_->calc_ctot_aq();

    const int noIter = solver_state_.noIterChargeBalance_;

    logger_->info("No global iter {}, final pH {:12.8e} and charge balance {:12.8e}",
                    noIter,
                    Vchem_->calc_pH(),
                    f_next);

    return has_converged;
}

void GCSolver::newton_ah_log(const RealVec& x, RealVec& Fchem, double** jacobi)
{
    ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;
	ChemTable &SM_mineral = *Vchem_->ICS_->SM_mineral_;

    const double F_div_RT = Vchem_->ICS_->CP_.F_div_RT_;
    const double kappa_s = calc_kappa();

    // Calculate E=exp(F*psi0/RT)
    double E{};
    if(Vchem_->surface_flag_)
    {
        const auto logE = static_cast<double>(x[Vchem_->size_mass_]);
        Vchem_->log_a_[Vchem_->pos_exp_] = logE;
        Vchem_->log_m_[Vchem_->pos_exp_] = logE;
        Vchem_->psi_ = logE * NumericalConstants::LNTEN / F_div_RT;
        E = POW10(logE);
    }

    // Construct basis vector
    for (int i = 0; i < Vchem_->size_mass_; ++i)
    {
        const int pos = Vchem_->pos_mass_[i];

        Vchem_->log_m_[pos] = static_cast<double>(x[i]);
        Vchem_->log_a_[pos] = Vchem_->log_m_[pos]+Vchem_->log_g_[pos];

        if (SM_basis.type_[pos] == GeochemicalComponentType::ION_EXCHANGE)
        {
            Vchem_->log_a_[pos] += log10(Vchem_->ctot_[pos]);
        }
    }

    for (int i = 0; i < Vchem_->size_gas_; ++i)
    {
        Vchem_->ICS_->SM_mineral_->log_af_[Vchem_->pos_gas_[i]] = static_cast<double>(x[Vchem_->size_mass_ + Vchem_->surface_flag_ + i]);
    }

    Vchem_->calc_rock_spec_conc();

    for(int i=0; i < Vchem_->size_; ++i)
    {
        mpow_[i] = POW10(Vchem_->log_m_[i]);  // m = 10^log_m_
    }

    // Calculate concentrations of complexes
    Vchem_->calc_complex_sparse(SM_all);
	// ... and SI indexes
    Vchem_->calc_complex_sparse(SM_mineral);
    for(int i=0; i < SM_all.noRows_; ++i)
    {
        smpow_[i] = POW10(SM_all.log_m_[i]);  // n = 10^log_SM_all_m
    }

    // Amount dissolved of rock-buffered minerals
    for(int i=0; i < Vchem_->size_rock_; ++i)
    {
        calc_F_spec(i, Vchem_->pos_rock_[i], mpow_, smpow_, delta_calc_);

        if(Vchem_->size_sup_min_>0)
        {
            calc_F_sup_min(i, Vchem_->pos_rock_[i], delta_calc_);
        }

        if (Vchem_->surface_options_.CALCULATE_DIFFUSE_LAYER)
        {
            // 18/7-22: Removed frac_DL_-factor.
            delta_calc_[i+1] -= Vchem_->ctot_dl_excess_[Vchem_->pos_rock_[i]];
        }
    }

    // Calculate F_j = c_tot_j - m_j - sum_i^N_rb M_ij delta_i-sum_i^Nx SM_i,j x_i.
    Fchem[0] = 0.0;
    for(int j=0; j < Vchem_->size_mass_ ; ++j)
    {
        // Basis & secondary species
        calc_F_spec(j, Vchem_->pos_mass_[j], mpow_, smpow_, Fchem);

        if (Vchem_->surface_options_.CALCULATE_DIFFUSE_LAYER)
        {
            // 18/7-22: Removed frac_DL_-factor.
            Fchem[j + 1] -= Vchem_->ctot_dl_excess_[Vchem_->pos_mass_[j]];
        }

        // ---> Add amounts lost or gained due to precipitation/dissolution

        if(Vchem_->size_sup_min_ > 0)
        {
            calc_F_sup_min(j, Vchem_->pos_mass_[j], Fchem);
        }
        if(Vchem_->size_rock_  > 0)
        {
            add_F_rock_min(j, Fchem, delta_calc_);
        }

        // ---> Add amount in a gas phase
        if (Vchem_->size_gase_phase_ > 0)
        {
            calc_F_gas_phase(j, Vchem_->pos_mass_[j], Fchem);
        }
    }

    if (Vchem_->surface_flag_)
    {
        // TODO: How to extend when using the Mean Potential Approximation?
        //       [a.k.a. "Donnan"] Do we need to solve a different equation?
        const auto& surf_options = Vchem_->surface_options_;
        if (surf_options.CALCULATE_DIFFUSE_LAYER && surf_options.INCLUDE_SURFACE_EXCESSES_IN_SURFACE_EQUATION)
        {
            calc_F_diffuse_layer_cbal_simple(Fchem, kappa_s, E, mpow_, smpow_);
        }
        else if(surf_options.GRAHAME_EQ_VERSION == GrahameEquation::GENERAL)
        {
            calc_F_grahame(Fchem, kappa_s, E, mpow_, smpow_);
        }
        else if (surf_options.GRAHAME_EQ_VERSION == GrahameEquation::SYMMETRICAL_ELECTROLYTE)
        {
            calc_F_grahame_symmetrical(Fchem, kappa_s, E, mpow_, smpow_);
        }
        else
        {
            assert(surf_options.GRAHAME_EQ_VERSION == GrahameEquation::ASSUME_CHARGE_BALANCE);
            calc_F_grahame_simple(Fchem, kappa_s, E, mpow_, smpow_);
        }
    }

    for (int j = 0; j < Vchem_->size_gas_; ++j)
    {
        calc_F_gas(j, Fchem);
    }

    if (!options_.CALC_JACOBI_NUM_ || Vchem_->ICS_->PRINT_DEBUG_CHEM_)
    {
        // Intermediate step for the rock part of the Jacobian
        if (Vchem_->size_rock_ > 0)
        {
            for (int j = 0; j < Vchem_->size_rock_; ++j)
            {
                const int pos_rock_j = Vchem_->pos_rock_[j];

                if (Vchem_->surface_flag_)
                {
                    calc_jacobi_spec(j, pos_rock_j, Vchem_->size_mass_, Vchem_->pos_exp_, Jacobi_dum_, mpow_, smpow_);
                }

                for (int q = 0; q < Vchem_->size_mass_; ++q)
                {
                    calc_jacobi_spec(j, pos_rock_j, q, Vchem_->pos_mass_[q], Jacobi_dum_, mpow_, smpow_);
                    calc_jacobi_sup_min(j, pos_rock_j, q, Vchem_->pos_mass_[q], Jacobi_dum_);
                }

                for (int q = 0; q < Vchem_->size_gas_; ++q)
                {
                    calc_jacobi_spec(j,
                                     pos_rock_j,
                                     Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                                     Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                                     Jacobi_dum_,
                                     mpow_,
                                     smpow_);

                    calc_jacobi_sup_min(j,
                                        pos_rock_j,
                                        Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                                        Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                                        Jacobi_dum_);
                }
            }
        }

        // jacobi_j,q = dF_j/dm_q = -delta_j,q m_j - ( sum_i^Nx SM_i,q * SM_i,j * x_i )
        for (int j = 0; j < Vchem_->size_mass_; ++j)
        {
            const int pos_mass_j = Vchem_->pos_mass_[j];

            for (int q = 0; q < Vchem_->size_mass_; ++q)
            {
                calc_jacobi_spec(j, pos_mass_j, q, Vchem_->pos_mass_[q], jacobi, mpow_, smpow_);

                if (Vchem_->size_rock_ > 0)
                {
                    calc_jacobi_rock_min(j, q, jacobi, Jacobi_dum_);
                }

                if (Vchem_->size_sup_min_ > 0)
                {
                    calc_jacobi_sup_min(j, pos_mass_j, q, Vchem_->pos_mass_[q], jacobi);
                }
				if (Vchem_->size_gase_phase_ > 0)
                {
                    calc_jacobi_gas_phase(j, pos_mass_j, q, Vchem_->pos_mass_[q], jacobi);
                }
            }

            if (Vchem_->surface_flag_)
            {
                calc_jacobi_spec(j,
                                 pos_mass_j,
                                 Vchem_->size_mass_,
                                 Vchem_->pos_exp_,
                                 jacobi,
                                 mpow_,
                                 smpow_);

                if (Vchem_->size_rock_ > 0)
                {
                    calc_jacobi_rock_min(j, Vchem_->size_mass_, jacobi, Jacobi_dum_);
                }

                const auto& surf_options = Vchem_->surface_options_;

                if (surf_options.CALCULATE_DIFFUSE_LAYER && surf_options.INCLUDE_SURFACE_EXCESSES_IN_SURFACE_EQUATION)
                {
                    calc_jacobi_diffuse_layer_cbal_simple_wrt_mq(j, pos_mass_j, jacobi, kappa_s, E, mpow_, smpow_);
                }
                else if(surf_options.GRAHAME_EQ_VERSION == GrahameEquation::GENERAL)
                {
                    calc_jacobi_grahame_wrt_mq(j, pos_mass_j, jacobi, kappa_s, E, mpow_, smpow_);
                }
                else if (surf_options.GRAHAME_EQ_VERSION == GrahameEquation::SYMMETRICAL_ELECTROLYTE)
                {
                    calc_jacobi_grahame_symmetrical_wrt_mq(j, pos_mass_j, jacobi, kappa_s, E, mpow_, smpow_);
                }
                else
                {
                    assert(surf_options.GRAHAME_EQ_VERSION == GrahameEquation::ASSUME_CHARGE_BALANCE);
                    calc_jacobi_grahame_wrt_mq_simple(j, pos_mass_j, jacobi, kappa_s, E, mpow_, smpow_);
                }
            }

            for (int q = 0; q < Vchem_->size_gas_; ++q)
            {
                calc_jacobi_spec(j,
                                 pos_mass_j,
                                 Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                                 Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                                 jacobi,
                                 mpow_,
                                 smpow_);

                calc_jacobi_rock_min(j,
                                     Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                                     jacobi,
                                     Jacobi_dum_);

                if (Vchem_->size_sup_min_ > 0)
                {
                    calc_jacobi_sup_min(j,
                                        pos_mass_j,
                                        Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                                        Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                                        jacobi);
                }
            }
        }

        if (Vchem_->surface_flag_)
        {
            const auto& surf_options = Vchem_->surface_options_;

            if (surf_options.CALCULATE_DIFFUSE_LAYER && surf_options.INCLUDE_SURFACE_EXCESSES_IN_SURFACE_EQUATION)
            {
                calc_jacobi_diffuse_layer_cbal_simple_wrt_E0(jacobi, kappa_s, E, mpow_, smpow_);

                for (int q = 0; q < Vchem_->size_gas_; ++q)
                {
                    calc_jacobi_diffuse_layer_cbal_simple_wrt_mq
                        (
                            Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                            Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                            jacobi,
                            kappa_s,
                            E,
                            mpow_,
                            smpow_
                        );
                }

            }
            else if (surf_options.GRAHAME_EQ_VERSION == GrahameEquation::GENERAL)
            {
                calc_jacobi_grahame_wrt_E0(jacobi, kappa_s, E, mpow_, smpow_);

                for (int q = 0; q < Vchem_->size_gas_; ++q)
                {
                    calc_jacobi_grahame_wrt_mq
                        (
                            Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                            Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                            jacobi,
                            kappa_s,
                            E,
                            mpow_,
                            smpow_
                        );
                }
            }
            else if (surf_options.GRAHAME_EQ_VERSION == GrahameEquation::SYMMETRICAL_ELECTROLYTE)
            {
                calc_jacobi_grahame_symmetrical_wrt_E0(jacobi, kappa_s, E, mpow_, smpow_);

                for (int q = 0; q < Vchem_->size_gas_; ++q)
                {
                    calc_jacobi_grahame_symmetrical_wrt_mq
                        (
                            Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                            Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                            jacobi,
                            kappa_s,
                            E,
                            mpow_,
                            smpow_
                        );
                }
            }
            else
            {
                assert(surf_options.GRAHAME_EQ_VERSION == GrahameEquation::ASSUME_CHARGE_BALANCE);
                calc_jacobi_grahame_wrt_E0_simple(jacobi, kappa_s, E, mpow_, smpow_);

                for (int q = 0; q < Vchem_->size_gas_; ++q)
                {
                    calc_jacobi_grahame_wrt_mq_simple(
                        Vchem_->size_mass_ + Vchem_->surface_flag_ + q,
                        Vchem_->pos_rock_[Vchem_->pos_rel_gas_[q]],
                        jacobi,
                        kappa_s,
                        E,
                        mpow_,
                        smpow_);
                }
            }
        } // end surface_flag_

        for (int j = 0; j < Vchem_->size_gas_; ++j)
        {
            for (int q = 0; q < Vchem_->size_mass_ + Vchem_->surface_flag_ + Vchem_->size_gas_; ++q)
            {
                calc_jacobi_gas(j, q, jacobi, Jacobi_dum_);
            }
        }
    }
}

void GCSolver::mass_balance()
{
    // ============================ Initialization =============================
    // Note: The initial values for ion exchange are really bad ...
    for(int i=0; i < Vchem_->size_mass_;++i)
    {
        const double log_mi = Vchem_->log_m_[Vchem_->pos_mass_[i]];
        x0_[i] = log_mi;
        x_init_[i] = log_mi;
    }

    const double E_init = (Vchem_->pos_exp_ > -1) ? Vchem_->log_a_[Vchem_->pos_exp_]: 0.0;
    x0_[Vchem_->size_mass_] = E_init;
    x_init_[Vchem_->size_mass_] = E_init;

    int dim = Vchem_->size_mass_ + Vchem_->surface_flag_;
    for (int i = 0; i < Vchem_->size_gas_; ++i)
    {
        const auto& log_af = Vchem_->ICS_->SM_mineral_->log_af_[Vchem_->pos_gas_[i]];
        x0_[dim + i] = log_af;
        x_init_[dim + i] = log_af;
    }
    dim += Vchem_->size_gas_;

    // ========================= Non-linear iterations =========================

    solver_state_.reset_mbal();

    double criterion;
    double delta;

    do
    {
        criterion = 0.0;

        if(options_.CALC_JACOBI_NUM_)
        {
            calc_analytical_and_numerical_jacobian(x_init_); // Writes file in run directory
        }
        else
        {
            newton_ah_log(x_init_, Fchem_, Jacobi_);
        }

        for (int i = 0; i < dim; ++i)
        {
            criterion += fabs(Fchem_[i + 1]);
        }

        if (Vchem_->ICS_->PRINT_DEBUG_CHEM_ == DebugInfoLevel::MAX)
        {
            write_mass_balance_debug_info_to_screen(DebugInfoTypeMBAL::BEFORE);
        }

        // If the residual equations are already satisfied, we stop here!
        if(criterion < options_.MBAL_CONV_CRITERION_)
        {
            break;
        }

        // Otherwise, we attempt to solve the linear system...
        int fale = ludcmp(Jacobi_, wksp_linalg_.data(), dim, chem_indx_.data(), Fchem_[0]);

        if(fale)
        {
            if (Vchem_->ICS_->PRINT_DEBUG_CHEM_)
            {
                write_jacobi(x_init_, Fchem_, Jacobi_, dim, "jmabl_fail.out");

                Vchem_->write("mbal_out.out");
                writeCurrentBasVecToFile("MBAL_FAILED.out");

                Vchem_->ICS_->SM_all_->write("mbal_SM_all.out");
                Vchem_->ICS_->SM_mineral_->write("mball_SM_mineral.out");
                Vchem_->ICS_->SM_basis_->write("mbal_SM_Basis.out");

                write_mass_balance_debug_info_to_screen(DebugInfoTypeMBAL::FAILED_LINSOLV);
            }

            solver_state_.massBalanceStatus_ = ConvergenceStatus::NOT_CONVERGED;
            return;
        }

        if(!fale)
        {
            lubksb(Jacobi_, dim, chem_indx_.data(), Fchem_.data());
        }

        const double MBAL_MAX_STEP = options_.MBAL_MAX_STEP_;
        for(int i=0; i<Vchem_->size_mass_; ++i)
        {
            delta = Fchem_[i + 1];
            if(!fale) delta = restrictByMaximumAbsoluteValue(delta, MBAL_MAX_STEP);
            else if(delta < 0) delta = -MBAL_MAX_STEP;
            else delta = MBAL_MAX_STEP;

            x_init_[i] -= delta;
            if(x_init_[i] > NumericalConstants::CONCENTRATION_THRESHOLD)
            {
                x_init_[i] = x0_[i];
            }
        }

        if (Vchem_->surface_flag_)
        {
            const double MBAL_SURF_MAX_STEP = options_.MBAL_SURF_MAX_STEP_;

            delta = Fchem_[Vchem_->size_mass_ + 1];
            if(!fale) delta = restrictByMaximumAbsoluteValue(delta, MBAL_SURF_MAX_STEP);
            else if(delta < 0) delta = -MBAL_SURF_MAX_STEP;
            else delta = MBAL_SURF_MAX_STEP;

            x_init_[Vchem_->size_mass_] -= delta;

            const double psi_init = x_init_[Vchem_->size_mass_] * NumericalConstants::LNTEN / Vchem_->ICS_->CP_.F_div_RT_;
            if (fabs(psi_init) >  NumericalConstants::POTENTIAL_THRESHOLD)
            {
                x_init_[Vchem_->size_mass_] = E_init;
            }
        }

        for (int i = 0; i < Vchem_->size_gas_; ++i)
        {
            delta = Fchem_[Vchem_->size_mass_ + Vchem_->surface_flag_ + i + 1];
            if(!fale) delta = restrictByMaximumAbsoluteValue(delta, MBAL_MAX_STEP);
            else if(delta < 0) delta = -MBAL_MAX_STEP;
            else delta = MBAL_MAX_STEP;

            x_init_[Vchem_->size_mass_ + Vchem_->surface_flag_ + i] -= delta;
        }

        ++solver_state_.noIterMassBalance_;
        if(solver_state_.noIterMassBalance_ > options_.CHEM_MAXITER_MASS_)
        {
            logger_->warning("MASS BALANCE DID NOT CONVERGE, FINAL VALUE = {:g}.", criterion);
            logger_->warning("NO MBAL ITERATIONS = {:d}, NO pH ITERATIONS = {:d}.", solver_state_.noIterMassBalance_, solver_state_.noIterChargeBalance_);
            logger_->info("Io={:g}\tpH={:g}\tTemp={:g}",
                          Vchem_->Io_,
                          Vchem_->calc_pH(),
                          Vchem_->Temp_[gFluidPhase::WATER] - 273.15);

            if(Vchem_->ICS_->PRINT_DEBUG_CHEM_ == DebugInfoLevel::MAX)
            {
                write_mass_balance_debug_info_to_screen(DebugInfoTypeMBAL::NO_CONVERGENCE);  // NB: Calls Newton again
            }

            solver_state_.massBalanceStatus_ = ConvergenceStatus::NOT_CONVERGED;
            return;
        }

    } while(criterion > options_.MBAL_CONV_CRITERION_);

    logger_->info("MBAL: No iter {:d} and final value {:12.8e}.", solver_state_.noIterMassBalance_, criterion);
    logger_->info("MBAL: Surface charge: {:g}\t Surface pot: {:g}.", Vchem_->sch_, Vchem_->psi_);
}

/**
 * Note: Solves the mass balance equations for a fixed pH.
 *
 * @param [in] x Guess for -pH  (the pH multiplied with -1).
 * @return Charge-balance error when using the input guess for -pH.
 */
double GCSolver::f_pH_charge_balance(double x)
{
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;
    ChemTable& SM_all = *Vchem_->ICS_->SM_all_;

    // Update activity and concentration of H+ (based on current guess)
    Vchem_->log_a_[Vchem_->pos_pH_] = x;
    Vchem_->log_m_[Vchem_->pos_pH_] = x - Vchem_->log_g_[Vchem_->pos_pH_];

    if (Vchem_->size_mass_ + Vchem_->size_gas_ > 0)
    {
        mass_balance();

        if (solver_state_.massBalanceStatus_ == ConvergenceStatus::NOT_CONVERGED)
        {
            return 1.0;  // FIXME: Should do something else here....
        }
    }
    else
    {
        // Calculate the concentration of the rock-buffered species and complexes
        Vchem_->calc_rock_spec_conc();

        for (int i = 0; i < Vchem_->size_; ++i)
        {
            mpow_[i] = POW10(Vchem_->log_m_[i]);  // m = 10^log_m_
        }

        // Update concentrations of complexes
        Vchem_->calc_complex_sparse(SM_all);

        for (int i = 0; i < SM_all.noRows_; ++i)
        {
            smpow_[i] = POW10(SM_all.log_m_[i]);  // n = 10^log_SM_all_m
        }
    }

    // Finally, we are ready to calculate the current charge-balance
    double charge_balance = 0.0;
    for (int i = 0; i < Vchem_->size_; ++i)
    {
        if (SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            charge_balance += mpow_[i] * SM_basis.charge_[i];

        }
    }
    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            charge_balance += smpow_[i] * SM_all.charge_[i];
        }
    }

    return charge_balance;
}

void GCSolver::add_F_rock_min(int j, RealVec& Fchem, const RealVec& delta_calc)
{
    ChemTable& SM_mineral = *Vchem_->ICS_->SM_mineral_;

    const int pos_j = Vchem_->pos_mass_[j];
    for (int i = 0; i < Vchem_->size_rock_; ++i)
    {
        const int pos_i = Vchem_->pos_buffer_[i];

        SM_mineral.delta_[pos_i] = 0.;

        for (int p = 0; p < Vchem_->size_rock_; ++p)
        {
            SM_mineral.delta_[pos_i] += Vchem_->sm_buf_[i][p] * delta_calc[p + 1];
        }

        Fchem[j + 1] -= SM_mineral.M_[pos_i][pos_j] * SM_mineral.delta_[pos_i];
    }

}

void GCSolver::calc_F_spec(int j, int pos, const RealVec& mpow, const RealVec& smpow, RealVec& Fchem) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    // Add basis species
    if(SM_basis.type_[pos] == GeochemicalComponentType::ION_EXCHANGE)
    {
        Fchem[j + 1] = Vchem_->ctot_[pos];
    }
    else
    {
        Fchem[j + 1] = Vchem_->ctot_[pos] - mpow[pos];
    }

    if (pos == Vchem_->pos_pH_)
    {
        Fchem[j + 1] -= Vchem_->WHTO_; // concentration of H2O
    }

    if (Vchem_->equilibrate_ && SM_basis.type_[pos] == GeochemicalComponentType::AQUEOUS_COMPLEX)
    {
        // Assume that c_tot is the bulk concentration and does not include
        // surface species. The mass balance is still valid for surface species.
        for (int i = 0; i < SM_all.noRows_; ++i)
        {
            if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
            {
                Fchem[j + 1] -= SM_all.M_[i][pos] * smpow[i];
            }
        }
    }
    else
    {
        for (int i = 0; i < SM_all.noRows_; ++i)
        {   // Add secondary species
            Fchem[j + 1] -= SM_all.M_[i][pos] * smpow[i];
        }
    }
}

/**
 * * NB: If rate equations are changed, remember to modify ChemBasVec::update_mineral_concentrations_kinetics() as well.
 *
 *   c_tot_calc = m_i + sum_j mu_ij n_i + e.q. minerals - (k_1+k_2*a_H)*(1-omega^m)^n
 *   F_j = c_tot-c_tot_calc */
void GCSolver::calc_F_sup_min(int j, int pos_j, RealVec& Fchem)
{
    const double aH = POW10(Vchem_->log_a_[Vchem_->pos_pH_]);

    for(int i=0; i < Vchem_->size_sup_min_; ++i)
    {
        const int posi = Vchem_->pos_sup_min_[i];

        // Calculate saturation index, log(IAP/K):
        double SI = 0.0;
        for(int k=0; k < Vchem_->ICS_->SM_mineral_->noColumns_; ++k)
        {
            SI += Vchem_->ICS_->SM_mineral_->M_[posi][k]*Vchem_->log_a_[k];
        }
        SI -= Vchem_->ICS_->SM_mineral_->logK_[posi];
        SI -= Vchem_->ICS_->SM_mineral_->log_af_[posi];

        Vchem_->ICS_->SM_mineral_->log_a_[posi] = SI;
        Vchem_->ICS_->SM_mineral_->log_m_[posi] = SI-Vchem_->ICS_->SM_mineral_->log_g_[posi];

        const double sgn = (SI < 0) ? 1. : -1.;
        const double k1_T = Vchem_->rate_laws_[i].k1_;
        const double k2_T = Vchem_->rate_laws_[i].k2_;
        const double m = Vchem_->rate_laws_[i].m_;
        const double n = Vchem_->rate_laws_[i].n_;
        const double n_acid = Vchem_->rate_laws_[i].n_acid_;
        const double stoichiometric_coefficient = Vchem_->ICS_->SM_mineral_->M_[posi][pos_j];

        // Note: Omega = 10^SI is the saturation state
        Fchem[j+1] += dt_*sgn*stoichiometric_coefficient*(k1_T + k2_T*pow(aH, n_acid))*pow(fabs(1.-POW10(SI*m)), n);
    }
}

/**
* We simply add amount of specie in gas phase
**/
void GCSolver::calc_F_gas_phase(int j, int pos_j, RealVec& Fchem)
{
    const double fact = 1e-3* Vchem_->gas_phase_frac_ / (PhysicalConstants::IdealGasConstant * Vchem_->Temp_[gFluidPhase::GAS]) * PhysicalConstants::standard_gas_pressure/Z_factor_gas_;

    for (int gas_i = 0; gas_i < Vchem_->size_gase_phase_; ++gas_i)
    {
        int pos_gas = Vchem_->pos_gas_phase_[gas_i];
        double SI = Vchem_->ICS_->SM_mineral_->log_a_[pos_gas];
        double fugacity = Vchem_->ICS_->SM_mineral_->fugacity_[pos_gas];
        Fchem[j + 1] -= Vchem_->ICS_->SM_mineral_->M_[pos_gas][pos_j] * POW10(SI) * fact / fugacity;
    }
}

void GCSolver::calc_F_gas(int j, RealVec& Fchem)
{
    const int phase = Vchem_->phase_[Vchem_->pos_rel_gas_[j]];
    const double Mw = Vchem_->ICS_->mol_weight_phase_[phase];
    const double frac = mass_phase_[gFluidPhase::WATER] / mass_phase_[phase];

    const int pos_gas = Vchem_->pos_gas_[j];

    double g_tmp = Vchem_->ctot_mineral_[pos_gas] + Vchem_->ICS_->SM_mineral_->delta_[pos_gas] * frac;
    g_tmp *= POW10(Vchem_->log_a_gas_[j]) * Mw;

    Fchem[Vchem_->size_mass_ + Vchem_->surface_flag_ + j + 1] = POW10(Vchem_->ICS_->SM_mineral_->log_af_[pos_gas]) - g_tmp;
}

void GCSolver::calc_F_grahame(RealVec& Fchem, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow)
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    calc_surface_charge(mpow, smpow);
    update_monovalent_divalent_charge_factor(E);  // 23/6-2022: Include this here also, so that we get n(1), n(2), etc. in the output..

    double g_tmp_1 = 0.0;
    for (int i = 0; i < Vchem_->size_; ++i)
    {
        if (SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            g_tmp_1 += mpow[i] * (pow(E, -SM_basis.charge_[i]) - 1.);
        }
    }

    double g_tmp_2 = 0.0;
    for (int i = 0; i < SM_all.noRows_; ++i) {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            g_tmp_2 += smpow[i] * (pow(E, -SM_all.charge_[i]) - 1.);
        }
    }

    Fchem[Vchem_->size_mass_ + 1] = kappa_s * (g_tmp_1 + g_tmp_2) - Vchem_->sch_ * Vchem_->sch_;
    Fchem[Vchem_->size_mass_ + 1] *= options_.MBAL_SCALE_SURF_;

}

void GCSolver::calc_F_grahame_simple(RealVec& Fchem, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow)
{
    calc_surface_charge(mpow, smpow);
    const double N = sqrt(update_monovalent_divalent_charge_factor(E));
    // Why not Fchem[Vchem_->size_mass_ + 1] = sqrt(kappa_s)* (sqrtE - 1.0 / sqrtE)*N; ????
    Fchem[Vchem_->size_mass_ + 1] = sqrt(4.0 * kappa_s) * sinh(0.5 * Vchem_->ICS_->CP_.F_div_RT_ * Vchem_->psi_) * N;
    Fchem[Vchem_->size_mass_ + 1] -= Vchem_->sch_;
    Fchem[Vchem_->size_mass_ + 1] *= options_.MBAL_SCALE_SURF_;
}

void GCSolver::calc_F_grahame_symmetrical(RealVec& Fchem, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow)
{
    // Note: This function assumes that there are only +1, -1, +2, -2 charged species
    calc_surface_charge(mpow, smpow);
    update_monovalent_divalent_charge_factor(E);

    const double sqrt_Io = sqrt(diffuse_layer_props_.calcIonicStrength());
    const double nu = Vchem_->surface_options_.VALENCE_OF_SYMMETRICAL_ELECTROLYTE;
    assert(nu == 1.0 || nu == 2.0);
    const double powE = std::pow(E, 0.5*nu); // before: sqrt(E);

    Fchem[Vchem_->size_mass_ + 1] = sqrt(kappa_s) * (powE - 1.0 / powE) * (sqrt_Io / nu);
    // Alternatively: Fchem[Vchem_->size_mass_ + 1] = (2.0/nu)*sqrt(kappa_s*I_o) * sinh(0.5 * nu * Vchem_->ICS_->CP_.F_div_RT_ * Vchem_->psi_);
    Fchem[Vchem_->size_mass_ + 1] -= Vchem_->sch_;
    Fchem[Vchem_->size_mass_ + 1] *= options_.MBAL_SCALE_SURF_;
}

void GCSolver::calc_jacobi_spec(int j, int pos, int q, int pos_q, double** jacobi, const RealVec& mpow, const RealVec& smpow) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_; // transformed matrix
    if (Vchem_->M_all_trans_.empty())
    {
        M_nu = SM_all.M_;
    }

    if(SM_basis.type_[pos] == GeochemicalComponentType::ION_EXCHANGE)
    {
        jacobi[j + 1][q + 1] = 0.0;
    }
    else
    {
        jacobi[j + 1][q + 1] = -Vchem_->beta_bas_inv_[pos][pos_q] * mpow[pos];
    }

    double g_tmp_q = 0.0;
    if (Vchem_->equilibrate_ && SM_basis.type_[pos] == GeochemicalComponentType::AQUEOUS_COMPLEX)
    {
        for (int i = 0; i < SM_all.noRows_; ++i)
        {
            if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
            {
                g_tmp_q += M_nu[i][pos_q] * SM_all.M_trans_[pos][i] * smpow[i];
            }
        }
    }
    else
    {
        for (int i = 0; i < SM_all.noRows_; ++i)
        {
            g_tmp_q += M_nu[i][pos_q] * SM_all.M_trans_[pos][i] * smpow[i];
        }
    }

    jacobi[j + 1][q + 1] -= g_tmp_q;
    jacobi[j + 1][q + 1] *= NumericalConstants::LNTEN;
}

void GCSolver::calc_jacobi_rock_min(int j, int q, double** jacobi, double** jac_dum) const
{
    const ChemTable& SM_mineral = *Vchem_->ICS_->SM_mineral_;
    const int pos = Vchem_->pos_mass_[j];

    double g_tmp_q = 0.0;
    for(int i=0; i < Vchem_->size_rock_; ++i)
    {
        const int pos_i = Vchem_->pos_buffer_[i];

        double g_tmp_p = 0.0;
        for(int p=0; p < Vchem_->size_rock_; ++p)
        {
            g_tmp_p += Vchem_->sm_buf_[i][p]*jac_dum[p+1][q+1];
        }
        g_tmp_q += SM_mineral.M_[pos_i][pos]*g_tmp_p;
    }
    jacobi[j+1][q+1] -= g_tmp_q;
}

void GCSolver::calc_jacobi_sup_min(int j, int pos, int q, int pos_q, double** jacobi) const
{
    const double aH = POW10(Vchem_->log_a_[Vchem_->pos_pH_]);

    double g_tmp_q = 0.0;

    for(int i=0; i < Vchem_->size_sup_min_; ++i)
    {
        const int pos_i  = Vchem_->pos_sup_min_[i];

        double dd;
        double dh;
        if(Vchem_->size_rock_ > 0)
        {
            dd = 0.0;
            for(int p=0; p < Vchem_->size_; ++p)
            {
                dd += Vchem_->ICS_->SM_mineral_->M_[pos_i][p]*Vchem_->beta_bas_inv_[p][pos_q];
            }

            dh = Vchem_->beta_bas_inv_[Vchem_->pos_pH_][pos_q];
        }
        else
        {
            dd = Vchem_->ICS_->SM_mineral_->M_[pos_i][pos_q];
            dh = 0.0;
        }

        const double k1 = Vchem_->rate_laws_[i].k1_;
        const double k2 = Vchem_->rate_laws_[i].k2_;
        const double m = Vchem_->rate_laws_[i].m_;
        const double n = Vchem_->rate_laws_[i].n_;
        const double n_acid = Vchem_->rate_laws_[i].n_acid_;

        const double SI = Vchem_->ICS_->SM_mineral_->log_a_[pos_i]; // log10(IAP/K)
        const double omega = POW10(SI); // IAP/K (saturation state)
        const double omega_pow = pow(omega, m);

        const double sgn = (SI < 0) ? 1.0 : -1.0;
        const double stoichiometric_coefficient = Vchem_->ICS_->SM_mineral_->M_[pos_i][pos];

        // Note: Uses that dOmega^m/dlog10(z) = ln(10)*z*dOmega^m/dz.
        const double d_dh = -k2*n_acid*pow(aH, n_acid)*pow(fabs(1.0-omega_pow), n)*dh;
        const double d_dd = (k1 + k2*pow(aH, n_acid)) * pow(fabs(1.0-omega_pow), n-1.0)*m*n*omega_pow*dd;
        const double derivative = dt_*sgn*stoichiometric_coefficient*(d_dh + d_dd);

        if(omega_pow < 1) g_tmp_q += derivative;
        else g_tmp_q -= derivative;
    }

    jacobi[j+1][q+1] -= g_tmp_q*NumericalConstants::LNTEN;
}

void GCSolver::calc_jacobi_gas_phase(int j, int pos, int q, int pos_q, double** jacobi) const
{
    double g_tmp_q = 0.0;
    const double fact = 1e-3* Vchem_->gas_phase_frac_ / (PhysicalConstants::IdealGasConstant * Vchem_->Temp_[gFluidPhase::GAS]) * PhysicalConstants::atmospheric_pressure/Z_factor_gas_;

    for(int gas_i=0; gas_i < Vchem_->size_gase_phase_; ++gas_i)
    {
        int pos_gas = Vchem_->pos_gas_phase_[gas_i];

        double dd;
        if(Vchem_->size_rock_ > 0)
        {
            dd = 0.0;
            for(int p=0; p < Vchem_->size_; ++p)
            {
                dd += Vchem_->ICS_->SM_mineral_->M_[pos_gas][p]*Vchem_->beta_bas_inv_[p][pos];
            }

        }
        else
        {
            dd = Vchem_->ICS_->SM_mineral_->M_[pos_gas][pos];
        }

        const double SI = Vchem_->ICS_->SM_mineral_->log_a_[pos_gas]; // log10(IAP/K)
        const double omega = POW10(SI); // IAP/K (saturation state)
        double fugacity = Vchem_->ICS_->SM_mineral_->fugacity_[pos_gas];

        g_tmp_q += fact*Vchem_->ICS_->SM_mineral_->M_[pos_gas][pos_q]*omega*dd / fugacity;

    }

    jacobi[j+1][q+1] -= g_tmp_q*NumericalConstants::LNTEN;
}
void GCSolver::calc_jacobi_gas(int j,  int q, double** jacobi, double** jac_dum) const
{
    const int pos_gj = Vchem_->pos_rel_gas_[j];
    const int phase = Vchem_->phase_[pos_gj];

    const double Mw = Vchem_->ICS_->mol_weight_phase_[phase];
    const double frac = mass_phase_[gFluidPhase::WATER] / mass_phase_[phase];

    const double g_tmp = POW10(Vchem_->log_a_gas_[j])*Mw*frac;

    double g_tmp_p = 0.0;
    for (int p = 0; p < Vchem_->size_rock_; ++p)
    {
        g_tmp_p += Vchem_->sm_buf_[pos_gj][p] * jac_dum[p + 1][q + 1]; // ln(10) included
    }

    double fact = 0.0;
    if (j + Vchem_->surface_flag_ + Vchem_->size_mass_ == q)
    {
        fact = POW10(Vchem_->ICS_->SM_mineral_->log_af_[Vchem_->pos_gas_[j]])*NumericalConstants::LNTEN;
    }

    jacobi[j+Vchem_->surface_flag_+Vchem_->size_mass_ + 1][q + 1] = fact - g_tmp * g_tmp_p;
}

void GCSolver::calc_jacobi_grahame_wrt_E0(double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const
{

    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if (Vchem_->M_all_trans_.empty())
    {
        M_nu = SM_all.M_;
    }

    // Start with contributions due to secondary species
    double g_tmp_1 = 0.0, g_tmp_2= 0.0, g_tmp_3=0.0;
    for(int i=0; i < SM_all.noRows_; ++i)
    {
        if(SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            const double Zi = SM_all.charge_[i];
            const double MiE0 = M_nu[i][Vchem_->pos_exp_];

            g_tmp_1 += Zi*smpow[i]*pow(E, -Zi);
            // (1/3-22): We have changed to -1 instead of -E here.
            g_tmp_2 += smpow[i]*MiE0*(pow(E, -Zi)-1.0); // <-- !
        }
        else if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            const double Zsci = SM_all.scharge_[i];
            const double MiE0 = M_nu[i][Vchem_->pos_exp_];
            g_tmp_3 += Zsci * smpow[i] * MiE0;
        }
    }


    //dn += SM_all->scharge_[i] * smpow[i] * M_nu[i][Vchem_->pos_exp_];
    static constexpr double F = PhysicalConstants::Faraday;
    const double SA = Vchem_->SA_;

    // Basis contribution
    jacobi[Vchem_->size_mass_+1][Vchem_->size_mass_+1] = kappa_s*(-g_tmp_1 + g_tmp_2) - 2.0*(F/SA)*Vchem_->sch_*g_tmp_3;

    g_tmp_1 = 0.0;
    for(int i=0; i<Vchem_->size_; ++i)
    {
        if(SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            g_tmp_1 += SM_basis.charge_[i]*mpow[i]*pow(E, -SM_basis.charge_[i]);
        }
    }

    jacobi[Vchem_->size_mass_+1][Vchem_->size_mass_+1] -= kappa_s*g_tmp_1;
    jacobi[Vchem_->size_mass_+1][Vchem_->size_mass_+1] *= NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;
}

void GCSolver::calc_jacobi_grahame_wrt_E0_simple(double** jacobi, double kappa_s, double E,
                                                 [[maybe_unused]] const RealVec& mpow, const RealVec& smpow) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if (Vchem_->M_all_trans_.empty())
    {
        M_nu = SM_all.M_;
    }

    const double sqrtE = sqrt(E);
    const double sinh2 = 0.5*(sqrtE - 1. / sqrtE); // sinh(beta*psi/2)
    const double cosh2 = 0.5*(sqrtE + 1. / sqrtE);
    const double sqrt_kappa = sqrt(kappa_s);

    double dn = 0.0, dn1 = 0.0, dn2 = 0.0, dn_2 = 0.0;
    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            const double dni = smpow[i] * M_nu[i][Vchem_->pos_exp_];
            const double Zi = SM_all.charge_[i];

            if (Zi > 1){
                dn2 += dni;
            }
            else if (Zi == 1){
                dn1 += dni;
            }
            else if (Zi < -1){
                dn_2 += dni;
            }
        }
        else
        if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            dn += SM_all.scharge_[i] * smpow[i] * M_nu[i][Vchem_->pos_exp_];
        }
    }

    static constexpr double F = PhysicalConstants::Faraday;
    const double SA = Vchem_->SA_;

    const double const_f = NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;
    const double N = sqrt(diffuse_layer_props_.calcIonicStrengthReplacementFactor(E));
    const double n2 = diffuse_layer_props_.n2_;
    const double n_2 = diffuse_layer_props_.n_2_;

    double der = sqrt_kappa*cosh2*N - (F/SA)*dn;
    der += sqrt_kappa*sinh2/N*(dn1 + dn2*(2.0 + 1.0/E) - n2/E + dn_2*E + n_2*E);
    jacobi[Vchem_->size_mass_ + 1][Vchem_->size_mass_ + 1] = const_f*der;
}

/** d(Grahame)/d log_10 m_q. */
void GCSolver::calc_jacobi_grahame_wrt_mq(int q, int pos_q, double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    static constexpr double F = PhysicalConstants::Faraday;
    const double SA = Vchem_->SA_;

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if (Vchem_->M_all_trans_.empty()){
        M_nu = SM_all.M_;
    }

    double dn = SM_basis.scharge_[pos_q] * mpow[pos_q]; // surface basis species never used in basis switching (never part of minerals)
    double g_tmp_1 = 0.0;
    for (int i = 0; i < SM_all.noRows_; ++i) // dsigma/dm_q
    {
        if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            dn += SM_all.scharge_[i] * smpow[i] * M_nu[i][pos_q];
        }
        else if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            g_tmp_1 += smpow[i] * M_nu[i][pos_q]*(pow(E, -SM_all.charge_[i]) - 1.0);
        }
    }
    double g_tmp_2 = 0.0;
    for (int i = 0; i < SM_basis.noRows_; ++i)
    {
        if (SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            g_tmp_2 += mpow[i] * Vchem_->beta_bas_inv_[i][pos_q]* (pow(E, -SM_basis.charge_[i]) - 1.0);
        }
    }

    jacobi[Vchem_->size_mass_ + 1][q + 1] = kappa_s * g_tmp_1 + kappa_s * g_tmp_2 - 2.0 * (F / SA) * Vchem_->sch_ * dn;
    jacobi[Vchem_->size_mass_+1][q+1] *= NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;

}

/** d(Grahame)/d log_10 m_q. */
void GCSolver::calc_jacobi_grahame_wrt_mq_simple(int q, int pos_q, double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    double dn1 = 0.0;
    double dn2 = 0.0;
    double dn_2 = 0.0;
    if (SM_basis.type_[pos_q] == GeochemicalComponentType::AQUEOUS_COMPLEX)
    {
        if (Vchem_->size_rock_ > 0)
        {
            for (int j = 0; j < Vchem_->size_; ++j)
            {
                const double g_tmp_1 = Vchem_->beta_bas_inv_[j][pos_q] * mpow[j];

                if (SM_basis.charge_[j] > 1) dn2 += g_tmp_1;
                else if (SM_basis.charge_[j] == 1) dn1 += g_tmp_1;
                else if (SM_basis.charge_[j] < -1) dn_2 += g_tmp_1;
            }
        }
        else
        {
            const double g_tmp_1 = mpow[pos_q];

            if (SM_basis.charge_[pos_q] > 1) dn2 = g_tmp_1;
            else if (SM_basis.charge_[pos_q] == 1) dn1 = g_tmp_1;
            else if (SM_basis.charge_[pos_q] < -1) dn_2 = g_tmp_1;
        }

    }

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if(Vchem_->M_all_trans_.empty()){
        M_nu = SM_all.M_;
    }

    double dn = SM_basis.scharge_[pos_q] * mpow[pos_q];
    for (int i = 0; i<SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            if (SM_all.charge_[i] > 1) dn2 += M_nu[i][pos_q] * smpow[i];
            else if (SM_all.charge_[i] == 1) dn1 += M_nu[i][pos_q] * smpow[i];
            else if (SM_all.charge_[i] < -1) dn_2 += M_nu[i][pos_q] * smpow[i];
        }
        else if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            dn += SM_all.scharge_[i] * smpow[i] * M_nu[i][pos_q];
        }
    }

    static constexpr double F = PhysicalConstants::Faraday;
    const double SA = Vchem_->SA_;

    const double const_f = NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;
    const double N = sqrt(diffuse_layer_props_.calcIonicStrengthReplacementFactor(E));
    const double Es = sqrt(E);

    double der = (0.5/N)*sqrt(kappa_s)*(Es - 1.0 / Es)*(dn1 + (2.0 + 1.0/E)*dn2 + E*dn_2) - (F/SA)*dn;
    jacobi[Vchem_->size_mass_ + 1][q + 1] = const_f*der;
}

/** IMPORTANT: This function still assumes that there are only +1, -1, +2, -2 species. */
void GCSolver::calc_jacobi_grahame_symmetrical_wrt_mq(int q, int pos_q, double** jacobi, double kappa_s, double E, const RealVec& mpow, const RealVec& smpow) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    double dn1 = 0.0;
    double dn_1 = 0.0;
    double dn2 = 0.0;
    double dn_2 = 0.0;

    if (SM_basis.type_[pos_q] == GeochemicalComponentType::AQUEOUS_COMPLEX)
    {
        if (Vchem_->size_rock_ > 0)
        {
            for (int j = 0; j < Vchem_->size_; ++j)
            {
                const double g_tmp_1 = Vchem_->beta_bas_inv_[j][pos_q] * mpow[j];

                if (SM_basis.charge_[j] > 1) dn2 += g_tmp_1;
                else if (SM_basis.charge_[j] == 1) dn1 += g_tmp_1;
                else if (SM_basis.charge_[j] == -1) dn_1 += g_tmp_1;
                else if (SM_basis.charge_[j] < -1) dn_2 += g_tmp_1;
            }
        }
        else
        {
            const double g_tmp_1 = mpow[pos_q];

            if (SM_basis.charge_[pos_q] > 1) dn2 = g_tmp_1;
            else if (SM_basis.charge_[pos_q] == 1) dn1 = g_tmp_1;
            else if (SM_basis.charge_[pos_q] == -1) dn_1 = g_tmp_1;
            else if (SM_basis.charge_[pos_q] < -1) dn_2 = g_tmp_1;
        }

    }

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if(Vchem_->M_all_trans_.empty())
    {
        M_nu = SM_all.M_;
    }

    double dn = SM_basis.scharge_[pos_q] * mpow[pos_q];
    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            if (SM_all.charge_[i] > 1) dn2 += M_nu[i][pos_q] * smpow[i];
            else if (SM_all.charge_[i] == 1) dn1 += M_nu[i][pos_q] * smpow[i];
            else if (SM_all.charge_[i] == -1) dn_1 += M_nu[i][pos_q] * smpow[i];
            else if (SM_all.charge_[i] < -1) dn_2 += M_nu[i][pos_q] * smpow[i];
        }
        else if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            dn += SM_all.scharge_[i] * smpow[i] * M_nu[i][pos_q];
        }
    }

    static constexpr double F = PhysicalConstants::Faraday;
    const double SA = Vchem_->SA_;
    const double sqrt_Io = sqrt(diffuse_layer_props_.calcIonicStrength());
    const double Io_der = 0.5 * (dn1 + dn_1) + 2.0 * (dn2 + dn_2);

    const double nu = Vchem_->surface_options_.VALENCE_OF_SYMMETRICAL_ELECTROLYTE;
    assert(nu == 1.0 || nu == 2.0);
    const double powE = std::pow(E, 0.5*nu); // before: sqrt(E);

    const double sinh2 = 0.5*(powE - 1.0 / powE);
    jacobi[Vchem_->size_mass_ + 1][q + 1] = (sinh2/sqrt_Io/nu)*sqrt(kappa_s)*Io_der - (F/SA)*dn;
    jacobi[Vchem_->size_mass_ + 1][q + 1] *= NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;

}

void GCSolver::calc_surface_charge(const RealVec& mpow, const RealVec& smpow)
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    Vchem_->sch_ = 0.0;
    for (int i = 0; i < Vchem_->size_; ++i)
    {
        if (SM_basis.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            Vchem_->sch_ += SM_basis.scharge_[i] * mpow[i];
        }
    }
    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            Vchem_->sch_ += SM_all.scharge_[i] * smpow[i];
        }
    }
    Vchem_->sch_ *= PhysicalConstants::Faraday / Vchem_->SA_;
}

/** @return The residual of the Grahame equation. Only works if called after mbal routine. */
double GCSolver::res_grahame_simple()
{
    calc_surface_charge(mpow_, smpow_);

    const double E = POW10(Vchem_->log_m_[Vchem_->pos_exp_]);
    // Note (22/12-21): Before, the class variables n1_, n2_, etc. were NOT updated by this function.
    const double N = sqrt(update_monovalent_divalent_charge_factor(E));
    const double kappa_s = calc_kappa();

    return Vchem_->sch_ - sqrt(4.*kappa_s)*sinh(0.5*Vchem_->ICS_->CP_.F_div_RT_*Vchem_->psi_)*N;
}

/** Calculates diffuse layer concentration. At exit, Vchem->ctot_dl contains the diffuse layer concentrations. */
void GCSolver::calc_diffuse_layer_conc(double* c_dl) const
{
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;

    if (Vchem_->log_m_[Vchem_->pos_exp_] == 0.0)
    {   // psi == 0
        for (int i = 0; i < Vchem_->size_; ++i)
        {
            c_dl[i] = 0.0;
        }
    }
    else
    {
        // First, calculate surface excesses of complexes
        std::vector<double> cc_DL;
        cc_DL.resize(SM_all.noRows_);

        for (int j = 0; j < SM_all.noRows_; ++j)
        {
            cc_DL[j] = 0.0;
            if (SM_all.type_[j] == GeochemicalComponentType::AQUEOUS_COMPLEX)
            {
                if (SM_all.charge_[j] < -1)
                    cc_DL[j] = diffuse_layer_props_.fDL_2_ * smpow_[j];
                else if (SM_all.charge_[j] == -1)
                    cc_DL[j] = diffuse_layer_props_.fDL_1_ * smpow_[j];
                else if (SM_all.charge_[j] == 1)
                    cc_DL[j] = diffuse_layer_props_.fDL1_ * smpow_[j];
                else if (SM_all.charge_[j] > 1)
                    cc_DL[j] = diffuse_layer_props_.fDL2_ * smpow_[j];
            }
        }

        for (int i = 0; i < Vchem_->size_; ++i)
        {
            if (SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
            {
                c_dl[i] = 0.0;
                if (SM_basis.charge_[i] < -1)
                    c_dl[i] = diffuse_layer_props_.fDL_2_ * mpow_[i];
                else if (SM_basis.charge_[i] == -1)
                    c_dl[i] = diffuse_layer_props_.fDL_1_ * mpow_[i];
                else if (SM_basis.charge_[i] == 1)
                    c_dl[i] = diffuse_layer_props_.fDL1_ * mpow_[i];
                else if (SM_basis.charge_[i] > 1)
                    c_dl[i] = diffuse_layer_props_.fDL2_ * mpow_[i];

                // Calculate complex contribution
                if (SM_basis.charge_[i] != 0.0)
                {
                    double g_tmp = 0.0;
                    for (int j = 0; j < SM_all.noRows_; ++j) {
                        g_tmp += SM_all.M_[j][i] * cc_DL[j];
                    }
                    c_dl[i] += g_tmp;
                }
            }
        }
    }
}

double GCSolver::calc_kappa() const
{
    static constexpr double Rg = PhysicalConstants::IdealGasConstant;
    static constexpr double e0 = PhysicalConstants::VacuumPermittivity;
    const double eW = Vchem_->ICS_->CP_.ew_;
    const double T = Vchem_->Temp_[gFluidPhase::WATER];
    return 2000.0*Rg*e0*eW*T;
}

double GCSolver::update_monovalent_divalent_charge_factor(double E)
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    diffuse_layer_props_.setConcentrationsToZero();

    for(int i = 0; i < SM_basis.noRows_; ++i)
    {
        if(SM_basis.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            if(SM_basis.charge_[i] > 1)
            {
                diffuse_layer_props_.n2_ += mpow_[i];
            }
            else if(SM_basis.charge_[i] == 1)
            {
                diffuse_layer_props_.n1_ += mpow_[i];
            }
            else if(SM_basis.charge_[i] == -1)
            {
                diffuse_layer_props_.n_1_ += mpow_[i];
            }
            else if(SM_basis.charge_[i] < -1)
            {
                diffuse_layer_props_.n_2_ += mpow_[i];
            }
        }
    }

    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            if (SM_all.charge_[i] > 1){
                diffuse_layer_props_.n2_ += smpow_[i];
            }
            else if(SM_all.charge_[i] == 1){
                diffuse_layer_props_.n1_ += smpow_[i];
            }
            else if(SM_all.charge_[i] == -1){
                diffuse_layer_props_.n_1_ += smpow_[i];
            }
            else if(SM_all.charge_[i] < -1){
                diffuse_layer_props_.n_2_ += smpow_[i];
            }
        }
    }

    return diffuse_layer_props_.calcIonicStrengthReplacementFactor(E);
}

void GCSolver::integrate_diffuse_layer_for_mono_or_divalent_ions()
{
    if (Vchem_->log_m_[Vchem_->pos_exp_] == 0.0)  // psi == 0
    {
        diffuse_layer_props_.setSurfaceExcessesToZero();
    }
    else
    {
        const double E = POW10(Vchem_->log_m_[Vchem_->pos_exp_]);

        switch (Vchem_->surface_options_.INTEGRATION_METHOD)
        {
        case DiffuseLayerIntegrationMethod::MEAN_POTENTIAL_APPROXIMATION: // a.k.a. "Donnan"
        {
            integrate_diffuse_layer_for_mono_or_divalent_ions_using_mean_potential_approximation(E);
            break;
        }
        case DiffuseLayerIntegrationMethod::AKSEL_HIORTH:
        {
            integrate_diffuse_layer_for_mono_or_divalent_ions_using_aksel_hiorth_method(E);
            break;
        }
        case DiffuseLayerIntegrationMethod::ANALYTICAL:
        {
            integrate_diffuse_layer_for_mono_or_divalent_ions_using_analytical_expressions(E);
            break;
        }
        case DiffuseLayerIntegrationMethod::BORKOVEC_WESTALL:
        {
            integrate_diffuse_layer_for_mono_or_divalent_ions_using_borkovec_westall_method(E);
            break;
        }
        default:
        {
            throw ModelDoesNotExistException("You need to select an integration method for the diffuse layer!");
            break;
        }
        }  // end switch statement

        if (Vchem_->surface_options_.ONLY_COUNTER_IONS && E != 1.0)
        {
            if (E < 1.0) diffuse_layer_props_.setSurfaceExcessesOfNegativeIonsToZero();
            else diffuse_layer_props_.setSurfaceExcessesOfPositiveIonsToZero();
        }
    } // end else (integration)
}

void GCSolver::integrate_diffuse_layer_for_mono_or_divalent_ions_using_aksel_hiorth_method(double E)
{
    const double x_i = 1.0;
    const double x_f = E;

    const double kappa_s = calc_kappa();
    const double gi_prefactor = +0.5 * (Vchem_->SA_ / PhysicalConstants::Faraday) * sqrt(kappa_s);

    diffuse_layer_props_.fDL_2_ = gi_prefactor * integrate_gi_ah(x_i, x_f, -2);
    diffuse_layer_props_.fDL_1_ = gi_prefactor * integrate_gi_ah(x_i, x_f, -1);
    diffuse_layer_props_.fDL1_ = gi_prefactor * integrate_gi_ah(x_i, x_f, 1);
    diffuse_layer_props_.fDL2_ = gi_prefactor * integrate_gi_ah(x_i, x_f, 2);
}

void GCSolver::integrate_diffuse_layer_for_mono_or_divalent_ions_using_borkovec_westall_method(double E)
{
    const double x_i = 1.0;
    const double x_f = E;

    const double kappa_s = calc_kappa();
    const double gi_prefactor = +0.5 * (Vchem_->SA_ / PhysicalConstants::Faraday) * sqrt(kappa_s);

    diffuse_layer_props_.fDL_2_ = gi_prefactor * integrate_gi_borkovec_westall(x_i, x_f, -2);
    diffuse_layer_props_.fDL_1_ = gi_prefactor * integrate_gi_borkovec_westall(x_i, x_f, -1);
    diffuse_layer_props_.fDL1_ = gi_prefactor * integrate_gi_borkovec_westall(x_i, x_f, 1);
    diffuse_layer_props_.fDL2_ = gi_prefactor * integrate_gi_borkovec_westall(x_i, x_f, 2);
}

void GCSolver::integrate_diffuse_layer_for_mono_or_divalent_ions_using_mean_potential_approximation(double E)
{
    // Assumes that we truncate the surface excess integral at x=d_DL ??
    const double gi_prefactor = 1.0e-6 * Vchem_->SA_ * Vchem_->d_DL_;
    const double E_inv = 1.0 / E;
    diffuse_layer_props_.fDL1_  = gi_prefactor*(E_inv - 1.0);
    diffuse_layer_props_.fDL2_  = gi_prefactor*(E_inv*E_inv - 1.0);
    diffuse_layer_props_.fDL_1_ = gi_prefactor*(E - 1.0);
    diffuse_layer_props_.fDL_2_ = gi_prefactor*(E*E - 1.0);
}

/** Analytical expressions using (essentially) the same method as in:
 *
 *   Oldham, K. B. (1975): "Composition of the diffuse double layer in seawater or other media containing ionic species of+ 2,+ 1,− 1, and− 2 charge types".
 *   Journal of Electroanalytical Chemistry and Interfacial Electrochemistry, 63(2), 139-156.
 */
void GCSolver::integrate_diffuse_layer_for_mono_or_divalent_ions_using_analytical_expressions(double E)
{
    diffuse_layer_props_.setSurfaceExcessesToZero();

    if(E == 1.0)
    {
        return;  // No surface excesses when psi == 0
    }

    static constexpr double ZERO_THRESHOLD = NumericalConstants::ZERO_CONCENTRATION_THRESHOLD;

    const double n1 = diffuse_layer_props_.n1_;
    const double n2 = diffuse_layer_props_.n2_;
    const double n_2 = diffuse_layer_props_.n_2_;

    const double nu = 1.0 - 1.0/E;
    const double kappa_s = calc_kappa();
    const double Io = n1 + 3.0*n2 + n_2;  // assumes charge-balance (eliminates -1 ions)
    const double n_1_from_chbal = n1 + 2.0 * (n2 - n_2); // ditto
    const double p = (n1 + 4.0*n2)/Io;
    const double q = n2/Io;
    const double gi_prefactor = +0.5 * (Vchem_->SA_ / PhysicalConstants::Faraday) * sqrt(kappa_s) / sqrt(Io);

    if (n2 > ZERO_THRESHOLD)  // q > 0
    {
        const double sqrt_q = sqrt(q);
        const double numer = p - 2.0*q*nu;
        const double denom = sqrt(p*p - 4.0*q); // TODO: When is this non-defined?
        const double acosh_f = acosh((2.0-p-nu*(p-2.0*q))/(1.0-nu)/denom);
        const double acosh_i = acosh((2.0-p)/denom);

        diffuse_layer_props_.fDL1_ = gi_prefactor * (1.0 / sqrt_q) * (acosh(numer / denom) - acosh(p / denom));
        diffuse_layer_props_.fDL_1_ = gi_prefactor * (acosh_f - acosh_i) / sqrt(q - p + 1.0);  // Note: (q - p + 1.0) = n(-2)/Io >= 0
    }
    else if (n_2 < ZERO_THRESHOLD)
    {
        // 1:1-electrolyte
        diffuse_layer_props_.fDL1_ = 2.0 * gi_prefactor * (sqrt(1.0 - nu) - 1.0);
        diffuse_layer_props_.fDL_1_ = 2.0 * gi_prefactor * (1.0 / sqrt(1.0 - nu) - 1.0);
    }
    else {
        // n2_ = 0, n_2_ > 0 (we always have n1_ > 0.0)
        const double fac = (1.0-p*nu)/(1.0-p);
        const double ln_argf = (sqrt(fac)+1)/(sqrt(fac)-1);
        const double ln_argi = (1.0 + sqrt(1.0-p))/(1.0 - sqrt(1.0-p));
        const double int_evalf = log(std::fabs(ln_argf));
        const double int_evali = log(std::fabs(ln_argi));

        // Can we simplify the above in terms, e.g., in terms of arccosh? (as for q > 0 case)
        diffuse_layer_props_.fDL1_ = gi_prefactor*2.0/p*(sqrt(1.0-p*nu) - 1.0);
        diffuse_layer_props_.fDL_1_ = gi_prefactor*(int_evalf-int_evali)/sqrt(1.0-p);
    }

    // RESOLVE: Multiply nu in the denominators below by an extra factor frac_DL? (think not..)
    if (n2 > ZERO_THRESHOLD)
    {
        //const double sigma_check = nu*sqrt(kappa_s*Io)*sqrt(q*nu*nu - p*nu + 1.0)/(1.0-nu);
        double g2_numer1 = Vchem_->sch_ * (1 - nu) / nu - sqrt(kappa_s * Io);
        g2_numer1 *= Vchem_->SA_ / PhysicalConstants::Faraday;
        diffuse_layer_props_.fDL2_ = 0.5 * (g2_numer1 - n1 * diffuse_layer_props_.fDL1_) / n2;
    }

    if (n_2 > ZERO_THRESHOLD)
    {
        double g_2_numer1 = Vchem_->sch_ / nu - sqrt(kappa_s * Io);
        g_2_numer1 *= Vchem_->SA_ / PhysicalConstants::Faraday;
        diffuse_layer_props_.fDL_2_ = 0.5 * (g_2_numer1 - n_1_from_chbal * diffuse_layer_props_.fDL_1_) / n_2;
    }
}


double GCSolver::integrate_gi_ah(double a, double b, double charge)
{
    // Set up callable to pass to the numerical integrator
    double (GCSolver::*fptr)(double, double) = &GCSolver::diffuse_layer_integrand_simple; // before: "Grahame_DL"
    auto f_integrand = [&fptr, this, charge](double E) {
      return (this->*fptr)(E, charge);
    };

    return numerical_integrator_.integrate(f_integrand, a, b);
}

double GCSolver::integrate_gi_borkovec_westall(double a, double b, double charge)
{
    double (GCSolver::*fptr)(double, double) = &GCSolver::diffuse_layer_integrand;
    auto f_integrand = [&fptr, this, charge](double E) {
      return (this->*fptr)(E, charge);
    };

    return numerical_integrator_.integrate(f_integrand, a, b);
}

/**
 * Essentially the integrand used in:
 *
 * Borkovec M, Westall J (1983) Solution of the Poisson-Boltzmann equation for surface
 * excesses of ions in the diffuse layer at the oxide-electrolyte interface.
 * J. Electroanal. Chem. 150:325-337.
 *
 * The main differences are:
 *      - We integrate wrt. E=1/X, instead of wrt. X.
 *      - We ASSUME no charged species of valence >= 3.
 */
double GCSolver::diffuse_layer_integrand(double E, double charge)
{
    if (E == 1)
    {
        return 0.0;
    }
    else
    {
        // Problem: Due to charge-balance errors occurring during Newton iterations,
        //          we may attempt to take the square root of something negative...
        const double denominator_pow2 = diffuse_layer_props_.calcSquareOfDiffuseLayerIntegrandDenominator(E);
        if (denominator_pow2 < 0)
        {
            throw IntegrandIsNotWellDefined("Attempting to take the square root of a negative number...");
        }
        const double sign = (E < 1) ? -1.0 : 1.0;
        const double denominator = sign*sqrt(denominator_pow2);
        if (charge > 1)
            return (1.0 / E / E - 1.0) / denominator;
        else if (charge == 1)
            return (1.0 / E - 1.0) / denominator;
        else if (charge == -1)
            return (E - 1.0) / denominator;
        else if (charge < -1)
            return (E*E - 1.0) / denominator;
        else
            return 0.0;
    }
}

/**
 * Simplified version of the diffuse layer integrand, specialized to the case with no
 * charged species of valence >= 3 and asssuming charge-balance to eliminate the
 * concentrations of species having charge -1.
 */
double GCSolver::diffuse_layer_integrand_simple(double E, double charge)
{
    if (E == 1)
    {
        return 0.0;
    }
    else
    {
        const double N = sqrt(diffuse_layer_props_.calcIonicStrengthReplacementFactor(E));
        const double sqrtE = sqrt(E);

        if (charge > 1)
            return -(E + 1.) / (E*E*sqrtE*N);
        else if (charge == 1)
            return -1. / (E*sqrtE*N);
        else if (charge == -1)
            return 1. / (sqrtE*N);
        else if (charge < -1)
            return (E + 1.) / (sqrtE*N);
    }
    return 0.0;
}

/** Note: A difference from the various impl. of the Grahame eq. is that we multiply both sides with S/F.*/
void GCSolver::calc_F_diffuse_layer_cbal_simple(RealVec& Fchem, [[maybe_unused]] double kappa_s,
                                                double E, const RealVec& mpow, const RealVec& smpow)
{
    calc_surface_charge(mpow, smpow);
    update_monovalent_divalent_charge_factor(E);

    const double SA =  Vchem_->SA_;
    static constexpr double F = PhysicalConstants::Faraday;

    Fchem[Vchem_->size_mass_ + 1]  = diffuse_layer_props_.calcTotalSurfaceExcessConcentration();
    Fchem[Vchem_->size_mass_ + 1] += (SA/F)*Vchem_->sch_;
    Fchem[Vchem_->size_mass_ + 1] *= options_.MBAL_SCALE_SURF_;
}

/** d(Grahame)/d log_10 m_q. */
void GCSolver::calc_jacobi_diffuse_layer_cbal_simple_wrt_mq(int q,
                                                            int pos_q,
                                                            double** jacobi,
                                                            [[maybe_unused]] double kappa_s,
                                                            [[maybe_unused]] double E,
                                                            const RealVec& mpow,
                                                            const RealVec& smpow) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;
    const ChemTable& SM_basis = *Vchem_->ICS_->SM_basis_;

    double dn1 = 0.0;
    double dn2 = 0.0;
    double dn_1 = 0.0;
    double dn_2 = 0.0;

    if (SM_basis.type_[pos_q] == GeochemicalComponentType::AQUEOUS_COMPLEX)
    {
        if (Vchem_->size_rock_ > 0)
        {
            for (int j = 0; j < Vchem_->size_; ++j)
            {
                const double g_tmp_1 = Vchem_->beta_bas_inv_[j][pos_q] * mpow[j];

                if(SM_basis.charge_[j] > 1) dn2 += g_tmp_1;
                else if(SM_basis.charge_[j] == 1) dn1 += g_tmp_1;
                else if(SM_basis.charge_[j] == -1) dn_1 += g_tmp_1;
                else if(SM_basis.charge_[j] < -1) dn_2 += g_tmp_1;
            }
        }
        else
        {
            const double g_tmp_1 = mpow[pos_q];
            if(SM_basis.charge_[pos_q] > 1) dn2 = g_tmp_1;
            else if(SM_basis.charge_[pos_q] == 1) dn1 = g_tmp_1;
            else if(SM_basis.charge_[pos_q] == -1) dn_1 = g_tmp_1;
            else if(SM_basis.charge_[pos_q] < -1) dn_2 = g_tmp_1;
        }
    }

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if (Vchem_->M_all_trans_.empty()) {
        M_nu = SM_all.M_;
    }

    double dn = SM_basis.scharge_[pos_q] * mpow[pos_q];

    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            if(SM_all.charge_[i] > 1) dn2 += M_nu[i][pos_q] * smpow[i];
            else if(SM_all.charge_[i] == 1) dn1 += M_nu[i][pos_q] * smpow[i];
            else if(SM_all.charge_[i] == -1) dn_1 += M_nu[i][pos_q] * smpow[i];
            else if(SM_all.charge_[i] < -1) dn_2 += M_nu[i][pos_q] * smpow[i];
        }
        else if(SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            dn += SM_all.scharge_[i] * smpow[i] * M_nu[i][pos_q];
        }
    }

    // The derivative below assumes dgi = 0
    const double const_f = NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;

    double der = diffuse_layer_props_.fDL1_ * dn1;
    der += 2.0 * diffuse_layer_props_.fDL2_ * dn2;
    der -= diffuse_layer_props_.fDL_1_ * dn_1;
    der -= 2.0 * diffuse_layer_props_.fDL_2_ * dn_2;
    // Until 18/7-22: Also multiplied "der" with frac_DL_...

    jacobi[Vchem_->size_mass_ + 1][q + 1] = const_f*(der + dn);
}

/** d(Grahame)/d log_10 E0. */
void GCSolver::calc_jacobi_diffuse_layer_cbal_simple_wrt_E0(double** jacobi,
                                                            [[maybe_unused]] double kappa_s,
                                                            [[maybe_unused]] double E,
                                                            [[maybe_unused]] const RealVec& mpow,
                                                            const RealVec& smpow) const
{
    // Note (3/1-22): Since derivatives of basis species wrt. to E0 vanish,
    //                we only need secondary species.
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if (Vchem_->M_all_trans_.empty()){
        M_nu = SM_all.M_;
    }

    double dn = 0.0;
    double dn1 = 0.0, dn2 = 0.0, dn_1 = 0.0, dn_2 = 0.0;

    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            const double dni = smpow[i] * M_nu[i][Vchem_->pos_exp_];

            if (SM_all.charge_[i] > 1) dn2 += dni;
            else if (SM_all.charge_[i] == 1) dn1 += dni;
            else if (SM_all.charge_[i] == -1) dn_1 += dni;
            else if (SM_all.charge_[i] < -1) dn_2 += dni;
        }
        else
        if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            dn += SM_all.scharge_[i] * smpow[i] * M_nu[i][Vchem_->pos_exp_];
        }
    }

    // Note: The derivative below assumes dgi = 0
    const double const_f = NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;

    double der = diffuse_layer_props_.fDL1_*dn1;
    der += 2.0*diffuse_layer_props_.fDL2_*dn2;
    der -= diffuse_layer_props_.fDL_1_*dn_1;
    der -= 2.0*diffuse_layer_props_.fDL_2_*dn_2;
    // Until 18/7-22: Multiplied "der" with frac_DL_ as well...
    // Before 28.05.2021: Subtracted (F/SA)*dn also from the expression below:
    jacobi[Vchem_->size_mass_ + 1][Vchem_->size_mass_ + 1] = const_f*(der + dn);
}

/** This function still assumes that there are only +1, -1, +2, -2 species. */
void GCSolver::calc_jacobi_grahame_symmetrical_wrt_E0(double** jacobi,
                                                      double kappa_s,
                                                      double E,
                                                      [[maybe_unused]] const RealVec& mpow,
                                                      const RealVec& smpow) const
{
    const ChemTable& SM_all = *Vchem_->ICS_->SM_all_;

    std::vector<std::vector<double > >& M_nu = Vchem_->M_all_trans_;
    if (Vchem_->M_all_trans_.empty())
    {
        M_nu = SM_all.M_;
    }

    double dn = 0.0;
    double dn1 = 0.0;
    double dn_1 = 0.0;
    double dn2 = 0.0;
    double dn_2 = 0.0;

    for (int i = 0; i < SM_all.noRows_; ++i)
    {
        if (SM_all.type_[i] == GeochemicalComponentType::AQUEOUS_COMPLEX)
        {
            const double dni = smpow[i] * M_nu[i][Vchem_->pos_exp_];
            const double Zi = SM_all.charge_[i];

            if (Zi > 1)             dn2 += dni;
            else if (Zi == 1)       dn1 += dni;
            else if (Zi == -1)      dn_1 += dni;
            else if (Zi < -1)       dn_2 += dni;
        }
        else
        if (SM_all.type_[i] == GeochemicalComponentType::SURFACE_COMPLEX)
        {
            dn += SM_all.scharge_[i] * smpow[i] * M_nu[i][Vchem_->pos_exp_];
        }
    }

    static constexpr double F = PhysicalConstants::Faraday;
    const double SA = Vchem_->SA_;

    const double nu = Vchem_->surface_options_.VALENCE_OF_SYMMETRICAL_ELECTROLYTE;
    assert(nu == 1.0 || nu == 2.0);
    const double powE = std::pow(E, 0.5*nu); // before: sqrt(E);

    const double sinh2 = 0.5*(powE - 1.0 / powE);  // Before: 0.5*(sqrtE - 1.0 / sqrtE). Also equal to: sinh(nu*beta*psi/2).
    const double cosh2 = 0.5*(powE + 1.0 / powE);  // before: 0.5*(sqrtE + 1.0 / sqrtE)
    const double sqrt_kappa = sqrt(kappa_s);
    const double sqrt_Io = sqrt(diffuse_layer_props_.calcIonicStrength());
    const double Io_der = 0.5 * (dn1 + dn_1) + 2.0 * (dn2 + dn_2);

    const double const_f = NumericalConstants::LNTEN*options_.MBAL_SCALE_SURF_;

    // NB: For the first term, we get an extra factor nu when taking the
    // derivative wrt. log10(E0), however it is cancelled by the extra factor
    // (1/nu) outside...
    double derivative = sqrt_kappa * cosh2 * sqrt_Io;
    // On the other hand, there is no such cancellation when we differentiate the ionic strength...
    derivative += sqrt_kappa * (sinh2 / sqrt_Io) * (Io_der / nu);
    derivative -= (F/SA)*dn;

    jacobi[Vchem_->size_mass_ + 1][Vchem_->size_mass_ + 1] = const_f*derivative;
}

void GCSolver::set_Jacobian_to_zero(){

    const int dim = Vchem_->dimension();

    for (int i = 0; i < dim; ++i){
        for (int j = 0; j < dim; ++j){
            Jacobi_[i + 1][j + 1] = 0.;
        }
    }

}

void GCSolver::writeCurrentBasVecToFile(const std::string& file_name)
{
    Vchem_->write_solution(file_name, &diffuse_layer_props_);
}

// Writes a table of logK values to file as a function of Temp at constant P
void GCSolver::write_logK_to_file_T(double T_start,
                                    double T_final,
                                    double dT,
                                    double P,
                                    const std::string& file_name,
                                    const int type)
{
    std::ofstream ofile(file_name, std::ios::out);

    ChemTable* CM;

    if (type == GeochemicalComponentType::MINERAL)
    {
        CM = Vchem_->ICS_->SM_mineral_.get();
    }
    else // has to be of type complex
    {
        CM = Vchem_->ICS_->SM_all_.get();
    }

    const auto N = static_cast<int>((T_final - T_start) / dT);
    double T = T_start;
    ofile << "T";
    for (int i = 0; i < CM->noRows_; ++i)
        if (CM->type_[i] == type)
            ofile << "\t" << CM->row_name_[i] ;
    ofile << "\n";

    for( int i = 0; i < N; ++i)
    {
        Vchem_->update_dG(*HKF_EOS_, T, P);
        Vchem_->update_logK(*HKF_EOS_, T);

        ofile << T;
        for (int j = 0; j < CM->noRows_; ++j)
        {
            if (CM->type_[j] == type) ofile << "\t" << CM->logK_[j];
        }
        ofile << "\n";
        T += dT;
    }
}

void GCSolver::write_jacobi(const RealVec& xinit, RealVec& Fchem, double** jc, int size, const char* fname)
{
    FILE* fp = my_fopen(fname, "w");

    for (int i = 0; i < Vchem_->ICS_->SM_all_->noRows_; ++i)
    {
        fmt::print(fp, "{:s}\t", Vchem_->ICS_->SM_all_->row_name_[i]);
    }
    fmt::print(fp, "\n");

    for (int i = 0; i < Vchem_->ICS_->SM_all_->noRows_; ++i)
    {
        fmt::print(fp, "{:g}\t", smpow_[i]);
    }
    fmt::print(fp, "\nlogKT:\n");
    for (int i = 0; i < Vchem_->ICS_->SM_all_->noRows_; ++i)
    {
        fmt::print(fp, "{:g}\t", Vchem_->ICS_->SM_all_->logK_[i]);
    }
    fmt::print(fp, "\n");

    for (int i = 0; i < Vchem_->ICS_->SM_mineral_->noRows_; ++i)
    {
        fmt::print(fp, "{:s}\t", Vchem_->ICS_->SM_mineral_->row_name_[i]);
    }
    fmt::print(fp, "\nlogKT:\n");
    for (int i = 0; i < Vchem_->ICS_->SM_mineral_->noRows_; ++i)
    {
        fmt::print(fp, "{:g}\t", Vchem_->ICS_->SM_mineral_->logK_[i]);
    }
    fmt::print(fp, "\n");

    if (Vchem_->pos_exp_ > -1)
    {
        fmt::print(fp, "E={:g}\n", Vchem_->log_a_[Vchem_->pos_exp_]);
    }

    for (int i = 0; i < size - 1; ++i)
    {
        fmt::print(fp, "{:s}\t", Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]]);
    }

    if (Vchem_->surface_flag_)
    {
        fmt::print(fp, "E\n");
    }
    else{
        fmt::print(fp, "{:s}\t", Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[size-1]]);
    }
    fmt::print(fp, "\nx_init:\n");
    for (int i = 0; i < size; ++i)
    {
        fmt::print(fp,"{:g}\t", xinit[i]);
    }
    fmt::print(fp, "\n");

    for (int i = 0; i < Vchem_->ICS_->SM_basis_->noRows_; ++i)
    {
        fmt::print(fp, "{:s}\t", Vchem_->ICS_->SM_basis_->row_name_[i]);
    }
    fmt::print(fp, "\nmpow:\n");
    for (int i = 0; i<size; ++i)
    {
        fmt::print(fp, "{:g}\t", mpow_[i]);
    }
    fmt::print(fp, "\nFchem:\n");
    for (int i = 0; i < size; ++i){
        fmt::print(fp,"{:g}\t", Fchem[i+1]);
    }
    fmt::print(fp, "\nJacobi:\n");
    for(int i=0; i < size; ++i)
    {
        for(int j=0; j < size; ++j)
        {
            fmt::print(fp,"{:4.8e}\t", jc[i+1][j+1]);
        }
        fmt::print(fp, "\n");
    }
    fclose(fp);
}

void GCSolver::calc_analytical_and_numerical_jacobian(RealVec& x)
{
    static constexpr double dx = 1e-7;

    set_Jacobian_to_zero();

    // Allocate temporary space...
    const int dim = Vchem_->dimension();

    std::vector<std::vector <double> > Jn;
    std::vector<std::vector <double> > Ja;
    std::vector<double> F, Fp;
    F.resize(dim);
    Fp.resize(dim);
    Jn.resize(dim);
    Ja.resize(dim);

    for (int i = 0; i < dim; ++i)
    {
        Ja[i].resize(dim);
        Jn[i].resize(dim);
    }

    // Analytical Jacobian
    newton_ah_log(x, Fchem_, Jacobi_);
    for (int i = 0; i < dim; ++i)
    {
        F[i] = Fchem_[i + 1];
    }

    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            Ja[i][j] = Jacobi_[i+1][j+1];
        }
    }

    // Numerical Jacobian
    for (int j = 0; j < dim; ++j)
    {
        x[j] += dx;
        newton_ah_log(x, Fchem_, Jacobi_);
        for (int i = 0; i < dim; ++i)
        {
            Jn[i][j] = (Fchem_[i + 1] - F[i]) / dx;
        }
        x[j] -= dx;
    }

    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            Jacobi_[i + 1][j + 1] = Jn[i][j];
        }
    }

    // Write output to file
    const std::string ofile = "j2_" + std::to_string(solver_state_.noIterMassBalance_) + ".out";
    std::ofstream off(ofile, std::ios::out);

    off << "Analytical Jacobian\n";
    for (int i = 0; i < dim; ++i)
    {
        if (i < Vchem_->size_mass_)
        {
            off << Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]] << "\t";
        }
        else if (Vchem_->surface_flag_ && i == Vchem_->size_mass_)
        {
            off << "E\t";
        }
        else{
            off << Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_gas_[i-Vchem_->size_mass_-Vchem_->surface_flag_]] << "\t";
        }
        for (int j = 0; j < dim; ++j)
        {
            off << Ja[i][j] << "\t";
        }
        off << "\n";
    }
    off << "\n";

    off << "Numerical Jacobian" << "\n";
    for (int i = 0; i < dim; ++i)
    {
        if (i < Vchem_->size_mass_)
        {
            off << Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]] << "\t";
        }
        else if(Vchem_->surface_flag_ && i == Vchem_->size_mass_)
        {
            off << "E\t";
        }
        else{
            off << Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_gas_[i - Vchem_->size_mass_ - Vchem_->surface_flag_]] << "\t";
        }
        for (int j = 0; j < dim; ++j)
        {
            off << Jn[i][j] << "\t";
        }
        off << "\n";
    }
}

double GCSolver::small_a_EOS(double Tc, double Pc)
{
    double fc = (PhysicalConstants::IdealGasConstant * Tc);
    return  0.457235 * fc * fc / Pc;
}

double GCSolver::small_b_EOS(double Tc, double Pc)
{
    double fc = (PhysicalConstants::IdealGasConstant * Tc);
    return  0.077796 * fc / Pc;
}

double GCSolver::alpha_EOS(double T, double omega, double Tc)
{
    double a=(1. + (0.37464 + 1.54226 * omega - 0.26992 * omega*omega) * (1 - sqrt(T / Tc)));
    return a * a;
}
double GCSolver::A_EOS(double P, double T, double omega, double Pc, double Tc)
{
    double f = (PhysicalConstants::IdealGasConstant *T);
    double a = small_a_EOS(Tc, Pc);
    return alpha_EOS(T, omega, Tc) * a * P / f / f;
}
double GCSolver::B_EOS(double P, double T, double Pc, double Tc)
{
    double f = (PhysicalConstants::IdealGasConstant * T);
    double b = small_b_EOS(Tc,Pc);
    return b * P / f;
}
// Redlich O. and Kwong J. (1949) On the thermodynamics of solutions.V.An equation of state.Fugacities of gaseous solutions.Chem.Rev. 44, 233–244.
// Peng D.-Y. and Robinson D. B. (1976) A new two-constant equation of state.Ind.Eng.Chem.Fundum. 15, 59–64.
// C. A. J. Appelo, D. L. Parkhurst, V. E. A. Post "Equations for calculating hydrogeochemical reactions of minerals and gases such as CO2 at high pressures and temperatures",  GCA 125, 2014
// Note: If gas mineral phases are entered by user - ideal gas law is used
void GCSolver::update_fugacity_and_mol_volume()
{

    for (int i = 0; i < Vchem_->size_gase_phase_; ++i)
    {
        double P = Vchem_->Pres_[gFluidPhase::GAS]; // NB this should be partial pressure if mixing formula not known
        double T = Vchem_->Temp_[gFluidPhase::GAS];
        int pos = Vchem_->pos_gas_phase_[i];
        double Tc = Vchem_->ICS_->SM_mineral_->Tc_[pos];
        double Pc = Vchem_->ICS_->SM_mineral_->Pc_[pos];
        Vchem_->ICS_->SM_mineral_->fugacity_[pos] = 1.; // default value
        Vchem_->ICS_->SM_mineral_->mol_volume_[pos] = PhysicalConstants::IdealGasConstant * T / P;// Ideal gas law m^3/mol
        if (Tc > 0 && Pc > 0) // if Peng Robinson EOS values are available
        {
            double omega = Vchem_->ICS_->SM_mineral_->omega_acc_[pos];
            double Ac = A_EOS(P, T, omega, Pc, Tc);
            double Bc = B_EOS(P, T, Pc, Tc);
            double upper_limit = 2;

            double param1 = 1.0;
            double param2 = -(1.0 - Bc);
            double param3 = (Ac - 2.0 * Bc - 3.0 * Bc * Bc);
            double param4 = -(Ac * Bc - Bc * Bc - Bc * Bc * Bc);
            auto CubSol = CubicRootSolver<double, ThrowOnError>::solve(param1, param2, param3, param4);
            double Zsol = CubSol[0];
            if (Zsol < std::abs(upper_limit) * .999)
            {
                Z_factor_gas_ = Zsol;
                // debug
                double fugacity = Zsol - 1 - std::log(Zsol - Bc) - Ac / (2.828 * Bc) * std::log((Zsol + 2.414 * Bc) / (Zsol - 0.414 * Bc));
                Vchem_->ICS_->SM_mineral_->fugacity_[pos] = std::exp(fugacity);
                Vchem_->ICS_->SM_mineral_->mol_volume_[pos] = Zsol * PhysicalConstants::IdealGasConstant * T / P;// m^3/mol
            } // no solution use ideal gas
            else
            {
                Z_factor_gas_ = Zsol;
                std::string warning_msg("Did not find a solution for Peng-Robinson EOS for gas phase: " + Vchem_->ICS_->SM_mineral_->row_name_[pos] + " using ideal gas law");
                logger_->warning(warning_msg.c_str());
            }
        }
    }
}
// stored in global vector alpha_a_gas_phase_tot_ and b_gas_phase_
void GCSolver::calculate_mixing_parameters_EOS()
{
    double binary_mixing_parameter = 0.2; // Whitson - mayby only valid for CO2?
    double kij = 0.;
    double T = Vchem_->Temp_[gFluidPhase::GAS];
    b_gas_phase_ = 0.;
    alpha_a_gas_phase_tot_ = 0.;
    for (int i = 0; i < Vchem_->size_gase_phase_; ++i)
    {
        int pos_i = Vchem_->pos_gas_phase_[i];
        double Tc_i = Vchem_->ICS_->SM_mineral_->Tc_[pos_i];
        double Pc_i = Vchem_->ICS_->SM_mineral_->Pc_[pos_i];
        double omega_i = Vchem_->ICS_->SM_mineral_->omega_acc_[pos_i];
        double y_i = Vchem_->mol_fraction_gas_phase_[i];

        b_gas_phase_ += Vchem_->mol_fraction_gas_phase_[i] *small_b_EOS(Tc_i,Pc_i);
        alpha_a_gas_phase_[i] = 0.;
        double a_a_alpha = small_a_EOS(Tc_i, Pc_i) * alpha_EOS(T, omega_i, Tc_i);
        for (int j = 0; j < Vchem_->size_gase_phase_; ++j)
        {

            if (i == j)
                kij = 0.;
            else
                kij = binary_mixing_parameter;
            int pos_j = Vchem_->pos_gas_phase_[j];
            double Tc_j = Vchem_->ICS_->SM_mineral_->Tc_[pos_j];
            double Pc_j = Vchem_->ICS_->SM_mineral_->Pc_[pos_j];
            double omega_j = Vchem_->ICS_->SM_mineral_->omega_acc_[pos_j];
            double fact= small_a_EOS(Tc_j, Pc_j) * alpha_EOS(T, omega_j, Tc_j) * a_a_alpha;
            double y_j = Vchem_->mol_fraction_gas_phase_[j];

            alpha_a_gas_phase_[i] += std::sqrt(fact) *  y_j*(1.-kij);
        }
        alpha_a_gas_phase_tot_ += alpha_a_gas_phase_[i]*y_i;
    }
}

// Redlich O. and Kwong J. (1949) On the thermodynamics of solutions.V.An equation of state.Fugacities of gaseous solutions.Chem.Rev. 44, 233–244.
// Peng D.-Y. and Robinson D. B. (1976) A new two-constant equation of state.Ind.Eng.Chem.Fundum. 15, 59–64.
// C. A. J. Appelo, D. L. Parkhurst, V. E. A. Post "Equations for calculating hydrogeochemical reactions of minerals and gases such as CO2 at high pressures and temperatures",  GCA 125, 2014
// Note: If gas mineral phases are entered by user - ideal gas law is used
void GCSolver::update_fugacity_and_mol_volume_mixtures()
{
    // maybe add a couple of extra iterations
    int MAXIT = 1;

    for (int extra_it = 0; extra_it < MAXIT; extra_it++)
    {
        Vchem_->Pres_[gFluidPhase::GAS] = 0.;
        for (int i = 0; i < Vchem_->size_gase_phase_; ++i) // gas phase pressure
        {
            int pos = Vchem_->pos_gas_phase_[i];
            double fi = POW10(Vchem_->ICS_->SM_mineral_->log_a_[pos]);
            Vchem_->gas_phase_pressure_[i] = fi / Vchem_->ICS_->SM_mineral_->fugacity_[pos] * PhysicalConstants::standard_gas_pressure;
            Vchem_->Pres_[gFluidPhase::GAS] += Vchem_->gas_phase_pressure_[i];
        }
        for (int i = 0; i < Vchem_->size_gase_phase_; ++i)
        {
            int pos = Vchem_->pos_gas_phase_[i];
            double fi = POW10(Vchem_->ICS_->SM_mineral_->log_a_[pos]) * PhysicalConstants::standard_gas_pressure; //fugacity
            Vchem_->mol_fraction_gas_phase_[i] = fi / (Vchem_->ICS_->SM_mineral_->fugacity_[pos] * Vchem_->Pres_[gFluidPhase::GAS]);
        }
        Vchem_->Pres_[gFluidPhase::OIL] = Vchem_->Pres_[gFluidPhase::WATER] = Vchem_->Pres_[gFluidPhase::GAS]; // ignore capillary pressure

        calculate_mixing_parameters_EOS();

        double P = Vchem_->Pres_[gFluidPhase::GAS];
        double T = Vchem_->Temp_[gFluidPhase::GAS];
        double f = 1. / (PhysicalConstants::IdealGasConstant * T);
        double upper_limit = 3;
        double lower_limit = .1;
        double Ac = alpha_a_gas_phase_tot_ * P * f * f;
        double Bc = b_gas_phase_ * P * f;
        auto EOS = [Ac, Bc](double Z)
            {
                return Z * Z * Z - (1 - Bc) * Z * Z + (Ac - 2 * Bc - 3 * Bc * Bc) * Z - (Ac * Bc - Bc * Bc - Bc * Bc * Bc);
            };



        int iter = 0;
        while (EOS(lower_limit) * EOS(upper_limit) > 0 && iter < 5)
        {
            upper_limit += 1;
            iter++;
        }

        double param1 = 1.0;
        double param2 = -(1.0 - Bc);
        double param3 = (Ac - 2.0 * Bc - 3.0 * Bc * Bc);
        double param4 = -(Ac * Bc - Bc * Bc - Bc * Bc * Bc);
        auto CubSol = CubicRootSolver<double, ThrowOnError>::solve(param1, param2, param3, param4);
        Z_factor_gas_ = CubSol[0];
        if (Z_factor_gas_ < std::abs(upper_limit) * .999)
        {
            for (int i = 0; i < Vchem_->size_gase_phase_; ++i)
            {
                int pos = Vchem_->pos_gas_phase_[i];
                double Tc = Vchem_->ICS_->SM_mineral_->Tc_[pos];
                double Pc = Vchem_->ICS_->SM_mineral_->Pc_[pos];
                double Br = small_b_EOS(Tc, Pc) / b_gas_phase_;

                double fugacity_coeff = Br * (Z_factor_gas_ - 1) - std::log(Z_factor_gas_ - Bc) - Ac / (2.828 * Bc) * (2 * alpha_a_gas_phase_[i] / alpha_a_gas_phase_tot_ - Br) * std::log((Z_factor_gas_ + 2.414 * Bc) / (Z_factor_gas_ - 0.414 * Bc));
                Vchem_->ICS_->SM_mineral_->fugacity_[pos] = std::exp(fugacity_coeff);

                //  if (Vchem_->ICS_->SM_mineral_->fugacity_[pos] > 1.)
                //      Vchem_->ICS_->SM_mineral_->fugacity_[pos] = .4;
                  // mole volume should not be calculated at current press T only at reference condition
                if (Vchem_->ICS_->SM_mineral_->model_[pos] == LogKModel::ANA)
                {
                    double fi = POW10(Vchem_->ICS_->SM_mineral_->log_a_[pos]) * PhysicalConstants::standard_gas_pressure; //fugacity
                    Vchem_->ICS_->SM_mineral_->mol_volume_[pos] = Z_factor_gas_ * Vchem_->ICS_->SM_mineral_->fugacity_[pos] * PhysicalConstants::IdealGasConstant * T / fi;// m^3/mol
                }
            }

        } // no solution use ideal gas
        else
        {
            std::string warning_msg("Did not find a solution for Peng-Robinson EOS for gas phase, using ideal gas law");

            logger_->warning(warning_msg.c_str());
            Z_factor_gas_ = 1.; // ideal gas solution
            for (int i = 0; i < Vchem_->size_gase_phase_; ++i)
            {
                int pos = Vchem_->pos_gas_phase_[i];
                [[maybe_unused]] double fi = POW10(Vchem_->ICS_->SM_mineral_->log_a_[pos]) * PhysicalConstants::standard_gas_pressure;
               // Vchem_->mol_fraction_gas_phase_[i] = fi / (Vchem_->ICS_->SM_mineral_->fugacity_[pos] * Vchem_->Pres_[gFluidPhase::GAS]);
            }
        }
    }

}



void GCSolver::write_mass_balance_debug_info_to_screen(DebugInfoTypeMBAL type)
{
    const int dim = Vchem_->dimension();

    switch(type){
    case DebugInfoTypeMBAL::BEFORE:
    {
        std::cout << "DEBUG INFO BEFORE MASS BALANCE CALCULATIONS..." << "\n";
        for (int i = 0; i < dim; ++i)
        {
            if (i < Vchem_->size_mass_){
                std::cout << Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]];
                std::cout << "\t" << Fchem_[i+1] << "\n";
            }
            else if (Vchem_->surface_flag_  && i == Vchem_->size_mass_){
                std::cout << "E\t" << Fchem_[i + 1] << "\n";
            }
            else{
                std::cout << Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_gas_[i - Vchem_->size_mass_ - Vchem_->surface_flag_]];
                std::cout << Fchem_[i + 1] << "\n";
            }
        }
        break;
    }
    case DebugInfoTypeMBAL::FAILED_LINSOLV:
    {
        fmt::print("Mass balance routine failed @ GC_Calls {:d} and ;iter_mbal_ {:d}, debug info in : mbal_out.out\n",
                   solver_state_.noCalls_,
                   solver_state_.noIterMassBalance_);
        fmt::print("Activity\tConcentration into\n");
        for (int i = 0; i < Vchem_->size_mass_; ++i)
        {
            fmt::print("{:s} = {:g}\t{:g}\n",
                       Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]],
                       POW10(x0_[i]),
                       Vchem_->ctot_[Vchem_->pos_mass_[i]]);
        }
        if (Vchem_->surface_flag_)
        {
            fmt::print("E = {:g}\t{:g}\n",
                       POW10(x0_[Vchem_->size_mass_]),
                       Vchem_->ctot_[Vchem_->pos_exp_]);
        }
        fmt::print("Mineral phases: ");
        for (int i = 0; i < Vchem_->size_sup_min_; ++i)
        {
            fmt::print("{:s} ({:g})\t",
                       Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_sup_min_[i]],
                       Vchem_->ctot_mineral_[Vchem_->pos_sup_min_[i]]);
        }
        fmt::print("Activity after\n");
        for (int i = 0; i < Vchem_->size_mass_; ++i)
        {
            fmt::print("{:s} = {:g}\n",
                       Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]],
                       POW10(x_init_[i]));
        }
        if (Vchem_->surface_flag_)
        {
            fmt::print("E = {:g}\n", POW10(x_init_[Vchem_->size_mass_]));
        }
        if (Vchem_->ICS_->PRINT_DEBUG_CHEM_ == DebugInfoLevel::ALOT)
        {
            Vchem_->write("mbal_fail_vchem.out");
            writeCurrentBasVecToFile("mbal_fail_solution.out");
        }

        break;
    }
    case DebugInfoTypeMBAL::NO_CONVERGENCE:
    {
        fmt::print("Concentration into\n");
        for (int i = 0; i < Vchem_->size_mass_; ++i)
        {
            fmt::print("{:s}: log_a={:g} C_tot={:g}\n",
                       Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]],
                       x0_[i],
                       Vchem_->ctot_[Vchem_->pos_mass_[i]]);
        }
        fmt::print("Concentration after\n");
        for (int i = 0; i < Vchem_->size_mass_; ++i)
        {
            fmt::print("{:s}: log_a={:g} C_tot={:g}\n",
                       Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]],
                       x_init_[i],
                       Vchem_->ctot_[Vchem_->pos_mass_[i]]);
        }
        for (int i = 0; i < Vchem_->size_rock_; ++i)
        {
            fmt::print("Rock Buffers: {:s} and basis {:s} \n",
                       Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_buffer_[i]],
                       Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_rock_[i]]);
        }

        for (int  i = 0; i < Vchem_->size_gas_; ++i)
        {
            const double mass_ratio = mass_phase_[gFluidPhase::WATER]
                                      / mass_phase_[gFluidPhase::OIL];

            fmt::print("{:s} delta {:g} CT {:g}\n",
                       Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_gas_[i]],
                       Vchem_->ICS_->SM_mineral_->delta_[Vchem_->pos_gas_[i]] * mass_ratio,
                       Vchem_->ctot_mineral_[Vchem_->pos_gas_[i]]);
        }
        writeCurrentBasVecToFile("SOLUTION_MBAL_MAXITER.out");
        Vchem_->write("VCHEM_MBAL_MAXITER.out");

        std::cout << "Cw " << mass_phase_[gFluidPhase::WATER];
        std::cout << " Co " << mass_phase_[gFluidPhase::OIL] << "\n";
        std::cout << "C_Co2_in " << Vchem_->ctot_mineral_[Vchem_->pos_gas_[0]] << "\n";

        newton_ah_log(x_init_, Fchem_, Jacobi_);

        for (int i = 0; i < Vchem_->size_rock_; ++i)
        {
            const double mass_ratio = mass_phase_[gFluidPhase::WATER]
                                      / mass_phase_[gFluidPhase::OIL];

            std::cout << Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_buffer_[i]] << " ";
            std::cout << Vchem_->ICS_->SM_mineral_->delta_[Vchem_->pos_buffer_[i]] * mass_ratio << "\n";
        }

        for (int i = 0; i < dim; ++i)
        {
            if (i < Vchem_->size_mass_)
            {
                std::cout << Vchem_->ICS_->SM_basis_->row_name_[Vchem_->pos_mass_[i]];
                std::cout << "\t" << Fchem_[i + 1] << "\n";
            }
            else if (Vchem_->surface_flag_ && i == Vchem_->size_mass_){
                std::cout << "E\t" << Fchem_[i + 1] << "\n";
            }
            else
            {
                std::cout << Vchem_->ICS_->SM_mineral_->row_name_[Vchem_->pos_gas_[i - Vchem_->size_mass_ - Vchem_->surface_flag_]];
                std::cout << Fchem_[i + 1] << "\n";
            }
        }
        break;
    }
    default:
        break;

    }
}

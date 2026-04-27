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
#include <opm/simulators/geochemistry/StandaloneSolvers.hpp>

#include <sstream>

#include <opm/simulators/geochemistry/IO/ErrorMsg.hpp>
#include <opm/simulators/geochemistry/IO/InputReader.hpp>
#include <opm/simulators/geochemistry/Thermo/ThermoTable.h>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

namespace {

std::vector<std::string> parseSpeciesBlock(const std::vector<std::string>& block_content)
{
    std::vector<std::string> species_names;
    for (const auto& line : block_content)
    {
        const auto tokens = split_string(line);
        if (!tokens.empty())
        {
            species_names.push_back(tokens.front());
        }
    }
    return species_names;
}

std::vector<double> parseNumericBlock(const std::vector<std::string>& block_content, const std::string& block_name)
{
    std::vector<double> values;
    for (const auto& line : block_content)
    {
        const auto tokens = split_string(line);
        for (const auto& token : tokens)
        {
            try
            {
                values.push_back(std::stod(token));
            }
            catch (const std::exception&)
            {
                error_and_exit("Could not parse value \"{:s}\" in block {:s}.", token, block_name);
            }
        }
    }
    return values;
}

std::vector<std::string> readSpeciesNames(InputReader& input_reader)
{
    const auto block_species = parseSpeciesBlock(input_reader.get_block_keyword_content("SPECIESLIST"));
    if (!block_species.empty())
    {
        return block_species;
    }

    const auto species_name = input_reader.get_simple_keyword_value("Species");
    if (species_name.empty())
    {
        error_and_exit("Need to specify at least one species using keyword Species or block SPECIESLIST.");
    }
    return {species_name};
}

std::vector<double> readScalarOrBlock(InputReader& input_reader,
                                      const std::string& simple_keyword,
                                      const std::string& block_keyword)
{
    const auto block_values = parseNumericBlock(input_reader.get_block_keyword_content(block_keyword), block_keyword);
    if (!block_values.empty())
    {
        return block_values;
    }

    try
    {
        return {std::stod(input_reader.get_simple_keyword_value(simple_keyword))};
    }
    catch (const std::exception&)
    {
        error_and_exit("Could not parse keyword {:s}.", simple_keyword);
    }

    return {};
}

void convertTemperaturesToCelsius(std::vector<double>& temperatures, const std::string& unit)
{
    const auto unit_upper = to_upper_case(unit);
    if (unit_upper == "C" || unit_upper == "CELSIUS")
    {
        return;
    }
    if (unit_upper == "K" || unit_upper == "KELVIN")
    {
        for (auto& value : temperatures)
        {
            value -= 273.15;
        }
        return;
    }

    error_and_exit("Unsupported TempUnit {:s}. Use C or K.", unit);
}

void convertPressuresToBar(std::vector<double>& pressures, const std::string& unit)
{
    const auto unit_upper = to_upper_case(unit);
    if (unit_upper == "BAR" || unit_upper == "BARS")
    {
        return;
    }
    if (unit_upper == "PA" || unit_upper == "PASCAL" || unit_upper == "PASCALS")
    {
        for (auto& value : pressures)
        {
            value *= 1.0e-5;
        }
        return;
    }

    error_and_exit("Unsupported PresUnit {:s}. Use Pa or bar.", unit);
}

bool parseBoolKeyword(const std::string& value, const std::string& keyword_name)
{
    try
    {
        return std::stoi(value) != 0;
    }
    catch (const std::exception&)
    {
        error_and_exit("Could not parse keyword {:s}.", keyword_name);
    }

    return false;
}

} // namespace

EquilibriumSolver::EquilibriumSolver()
    : CGeoChemIF()
{
}

/** A simple routine for doing equilibrium calculations. */
BasVecInfo EquilibriumSolver::solve(const std::string& case_name,
                                    std::istream& inputStream)
{
    std::map<std::string, std::string> simple_key_value_pairs = std::map<std::string, std::string>();
    simple_key_value_pairs["Temp"] = "25.0";
    simple_key_value_pairs["Pres"] = "1.0e5";
    simple_key_value_pairs["equilibrate"] = "0";
    simple_key_value_pairs["fDL"] = "1";
    simple_key_value_pairs["debug"] = "0";
    simple_key_value_pairs["serialize"] = "0";

    InputReader header(simple_key_value_pairs, false); // case-insensitive keywords
    header.read_key_value_pairs(inputStream);

    // NB: Must read "ADD_SPECIES", etc. & correctly set database also for equilibrium calculations...
    ptrInputReader_->read(inputStream);
    set_geochemical_database_modifications_from_user_input();

    geoChemPhases_.resetFromInputStream(inputStream);

    const auto& all_names_and_types = geoChemPhases_.getPhaseNamesAndTypes();

    // TO DO: Use ICS_full_ here as well? (for consistency + reuse code)
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const int max_size = geoChemPhases_.selectPhases(geoChemPhases_.getPhaseNamesAndTypes(), phases_to_be_used);

    auto ICS_ptr = InitChem::create_from_input_data
        (
            nullptr,
            db_changes_made_by_user_,
            phases_to_be_used,
            max_size,
            "EQ_SOLVER"
        );

    ICS_ptr->Temp_ = 273.15 + std::stod(header.get_simple_keyword_value("Temp"));
    ICS_ptr->Pres_ = std::stod(header.get_simple_keyword_value("Pres"));
    const auto equilibrate = static_cast<bool>(std::stoi(header.get_simple_keyword_value("equilibrate")));
    const int DEBUG_INFO = std::stoi(header.get_simple_keyword_value("debug"));
    const auto SERIALIZE = static_cast<bool>(std::stoi(header.get_simple_keyword_value("serialize")));

    ICS_ptr->PRINT_DEBUG_CHEM_ = DEBUG_INFO;
    if (DEBUG_INFO)
    {
        logger_.setLogLevel(LogLevel::debug);
    }

    const double frac_DL = std::stod(header.get_simple_keyword_value("fDL"));

    std::vector<double> Ct;
    for (int i = 0; i < ICS_ptr->size_aq_; ++i)
    {
        ICS_ptr->kcmap_[i] = i;
        Ct.push_back(ICS_ptr->c_vchem_[i]);
    }

    const int firstIonExchangerIndex = ICS_ptr->size_aq_;
    for (int i = 0; i < ICS_ptr->size_io_; ++i)
    {
        ICS_ptr->kcmap_[firstIonExchangerIndex + i] = firstIonExchangerIndex + i;
        ICS_ptr->c_vchem_[firstIonExchangerIndex + i] = ICS_ptr->c_io_[i];
        Ct.push_back(ICS_ptr->c_io_[i]);
    }

    const double specific_surface_area = ICS_ptr->SA_;

    const int firstSurfaceComplexIndex = ICS_ptr->size_aq_ + ICS_ptr->size_io_;
    for (int i = 0; i < ICS_ptr->size_surf_; ++i)
    {
        const auto surface_molar_conc = (1.0e18 / PhysicalConstants::Avogadro) * specific_surface_area * ICS_ptr->c_surf_[i];
        ICS_ptr->kcmap_[firstSurfaceComplexIndex + i] = firstSurfaceComplexIndex + i;
        ICS_ptr->c_vchem_[firstSurfaceComplexIndex + i] = surface_molar_conc;
        Ct.push_back(surface_molar_conc);
    }

    for (int i = 0; i < ICS_ptr->size_min_; ++i)
    {
        ICS_ptr->krmap_[i] = i;
    }

    std::vector<double> Cads(ICS_ptr->size_basis_, 0.0);
    std::vector<double> scharge(2);
    double pH;

    BasVecInfo info;
    equilibrate_solution(ICS_ptr.get(),
                         equilibrate,
                         case_name,
                         DEBUG_INFO,
                         Ct,
                         Cads,
                         pH,
                         scharge,
                         &info,
                         frac_DL);


    if (SERIALIZE)
    {
        // Store contents of BasVecInfo in a special file to be able to reconstruct it later (for testing)
        std::ofstream cereal_file(case_name + ".cereal", std::ios::binary);
        cereal::BinaryOutputArchive archive(cereal_file);
        archive(info);
    }

    std::cout << " Final pH=" << pH << " Final Surface Charge " << scharge[0] << "\n";
    return info;
}


OneDimensionalTransportSolver::OneDimensionalTransportSolver()
    : CGeoChemIF()
{
}

/** A simple 1D transport solver, currently without molecular diffusion
*  For testing purposes, optionally returns a vector with the effluent ion
*   concentrations (default: empty vector).
*/
std::vector<EffluentIonData> OneDimensionalTransportSolver::solve(const std::string& case_name,
                                                                  std::istream& inputStream)
{
    // Create instance of InputReader local to this function - used to read keywords only applicable to solve()...
    std::map<std::string, std::string> simple_key_value_pairs = std::map<std::string, std::string>();
    simple_key_value_pairs["Tf"] = "24.0"; // Total simulation time (hours)
    simple_key_value_pairs["VolRate"] = "0.1"; // Flow rate ml/min
    // Eventually replace these by "inject" keyword..?

    simple_key_value_pairs["Temp"] = "25.0";
    simple_key_value_pairs["Pres"] = "10e5";
    simple_key_value_pairs["NoBlocks"] = "10";
    // Change (18/8-22): By default, assume a porosity phi=0.25 (before: 1.0).
    simple_key_value_pairs["Porosity"] = "0.25"; // Constant porosity
    simple_key_value_pairs["Volume"] = "25.0"; // System pore volume (mL)
    simple_key_value_pairs["Imp"] = "1"; // Implicit or explicit integration?
    simple_key_value_pairs["Flush"] = "0.5"; // Volume flushed out of each cell

    simple_key_value_pairs["punch"] = "-1"; // writes out block
    simple_key_value_pairs["WriteNetCDF"] = "0";
    simple_key_value_pairs["serialize"] = "0";

    InputReader header(simple_key_value_pairs, false);
    header.read(inputStream); // NB: reads keywords before geochem

    // For the code below to work, we require that all keywords like "add_species", "chemtol", "debug", etc. are placed either...
    //
    //  a) before the "GEOCHEM" keyword, or
    //  b) after ALL the "phase keywords" (rate, iexchange, solution, etc.)
    ptrInputReader_->read(inputStream);

    const int debug = std::stoi(ptrInputReader_->get_simple_keyword_value("DEBUG"));
    const int interpolate = std::stoi(ptrInputReader_->get_simple_keyword_value("INTERPOLATE"));
    set_geochemical_database_modifications_from_user_input();

    geoChemPhases_.resetFromInputStream(inputStream);

    // Read input parameters & convert to internal units
    const double uf_rate = 1.0e-3 / 60.0;
    const double Tf_seconds = 3600.0 * std::stod(header.get_simple_keyword_value("Tf")); // sec (input in hours)
    const double InitialVolRate = uf_rate*std::stod(header.get_simple_keyword_value("VolRate")); // L/sec

    const double Temp = 273.15 + std::stod(header.get_simple_keyword_value("Temp")); // K (input in C)
    const double Pres = std::stod(header.get_simple_keyword_value("Pres"));  // Pa (same as input)
    const double Volume = 1e-3 * std::stod(header.get_simple_keyword_value("Volume")); // L (input in mL)
    const double Flush = std::stod(header.get_simple_keyword_value("Flush"));

    double porosity = std::stod(header.get_simple_keyword_value("Porosity"));
    if(porosity <= 0.0 || porosity > 1.0)
    {
        error_and_exit("The input porosity needs to be in the range (0, 1].");
    }
    const int NoBlocks = std::stoi(header.get_simple_keyword_value("NoBlocks"));
    const double BlockVolume = Volume / NoBlocks;  // Block volume assumed constant.

    const int punch = std::stoi(header.get_simple_keyword_value("punch"));
    const int TRY_TO_WRITE_NETCDF = std::stoi(header.get_simple_keyword_value("WriteNetCDF"));

    // 1: store BasVecInfo objects for last block, 2: also create .cereal archive file
    const int SERIALIZE = std::stoi(header.get_simple_keyword_value("serialize"));

    int Imp = std::stoi(header.get_simple_keyword_value("Imp"));
    if (Imp == 0 && Flush > 1)
    {
        std::cout << "Cannot use explicit integration scheme if more than one pore volume is flushed to a block.. \n";
        std::cout << "..continue with implicit integration scheme!\n";
        Imp = 1;
    }

    // Read when and what to inject ------->
    SimulationSchedule simulation_schedule;

    std::vector<std::string> unique_names_inject = ptrInputReader_->unique_names("INJECT");
    std::vector<std::string> unique_names_solution = geoChemPhases_.getAllNamesOfType(GeochemicalPhaseType::AQUEOUS_SOLUTION);

    // CHANGE (18/8-22):
    //  - We now require both the solution id and dt to be defined in the block.
    //  - Each INJECT keyword should be given a unique id, which is not used for
    //    any other purpose than to distinguish between them.
    double previous_flow_rate = InitialVolRate;

    for(const auto& inject_keyword_occurrence: unique_names_inject)
    {
        const auto& content_as_vec = ptrInputReader_->get_block_keyword_content(to_upper_case(inject_keyword_occurrence));
        if(content_as_vec.size() != 1)
        {
            error_and_exit("Cannot parse INJECT keyword, there should only be a single line holding the solution id and the time step");
        }
        assert(content_as_vec.size() == 1);

        std::string solution_name;
        std::istringstream iss(content_as_vec[0]);
        iss >> solution_name;
        const std::string full_solution_name = "SOLUTION " + solution_name;
        int solution_index = index_of_element_in_container(full_solution_name, unique_names_solution);
        if(solution_index == -1)
        {
            error_and_exit("The solution \"{:s}\" given in keyword INJECT does not exist.", solution_name);
        }

        double dt;
        auto ret_dt = parse_value_or_error<double>(iss, "Could not parse time step for solution \"" + solution_name + "\" in keyword INJECT.\n");
        if(!ret_dt) std::exit(-1);
        else dt = ret_dt.value();

        double flow_rate;
        auto ret_q = parse_value_or_error<double>(iss, "Missing flow rate from keyword INJECT, use previous one.\n");
        if(!ret_q && previous_flow_rate > 0) flow_rate = previous_flow_rate;
        else if(!ret_q) std::exit(-1);
        else
        {
            flow_rate = uf_rate*ret_q.value();
            previous_flow_rate = flow_rate;
        }

        simulation_schedule.addEpoch(3600.0 * dt, GeochemBoundaryCondition(solution_index, flow_rate));
    }

    if(!simulation_schedule.isValid())
    {
        // If no INJECT keyword is given, assume a single epoch from t=0 to Tf,
        // and that we are always injecting the solution with id "1".
        assert(geoChemPhases_.hasPhaseWithName(GeochemicalPhaseType::AQUEOUS_SOLUTION, "SOLUTION 1"));
        simulation_schedule.addEpoch(Tf_seconds, GeochemBoundaryCondition(1, InitialVolRate));
    }
    // <------- read when and what to inject

    // ------------- Initialize ---------------------
    const int nbl = NoBlocks;

    // NB: Must set interpolate flag before allocating memory for solver, splay tree, etc.
    create_ICS_full();
    ICS_full_->INTERPOLATE_ = interpolate;
    ICS_full_->Pres_ = Pres;
    ICS_full_->Temp_ = Temp;

    ICS_full_->PRINT_DEBUG_CHEM_ = debug;
    if (debug)
    {
        logger_.setLogLevel(LogLevel::debug);
    }

    ICS_full_->possibly_change_inconsistent_options();

    allocate_memory_for_solver_and_splay_tree_etc();
    modify_selected_solver_options(*GCS_);

    // Initialize by assuming the first solution is in the reservoir, and the remaining ones are to be injected (possibly)
    if (!ICS_sol_.empty()) throw MultipleInitializationException("CGeoChemIF::solve(): Have already allocated memory for ICS_sol_...");

    auto [reservoir_phases, names_of_injected_solutions] = geoChemPhases_.GetInjectReservoirSimple();
    AddPhases(reservoir_phases, "RESERVOIR");
    for (const auto& solution_name : names_of_injected_solutions)
    {
        AddAqueousSolution(solution_name);
    }

    // Make some modifications to the InitChem objects...
    // (the first one is assumed to be the initial solution)
    ICS_sol_[0]->charge_balance_ = ICS_full_->charge_balance_;
    ICS_sol_[0]->Temp_ = Temp;
    ICS_sol_[0]->Pres_ = Pres;
    for (std::size_t solutionIdx = 1; solutionIdx < ICS_sol_.size(); ++solutionIdx)
    {
        // NB: We always assume charge-balance for injected solutions.
        ICS_sol_[solutionIdx]->charge_balance_ = true;
        ICS_sol_[solutionIdx]->Temp_ = Temp;
        ICS_sol_[solutionIdx]->Pres_ = Pres;
    }

    // Declare a lot of local variables, only to be used by the 1D transport
    // solver...
    std::vector<std::vector<double>> Ctot;
    std::vector<std::vector<double>> Cads;
    std::vector<std::vector<double>> Caq;
    std::vector<std::vector<double>> Caq_old;
    std::vector<std::vector<double>> Cmin;
    std::vector<std::vector<double>> C_init;
    std::vector<std::vector<double>> scharge;
    std::vector<std::vector<double>> C_inj;
    std::vector<std::vector<double>> C_DL;

    std::vector<double> C_ads_init;
    std::vector<double> pH;
    std::vector<double> cmin_inj;
    std::vector<double> log_min_inj;
    std::vector<double> Ceff;

    // If NetCDF file is to be written, add individual species concentrations to the output.
    // We only store values at the latest time step, and currently we skip t=0.
    std::map<std::string, std::vector<double>> block_concentrations_for_individual_species;
    std::map<std::string, std::vector<double>> block_precipitation_for_minerals; // NOTE: cumulative amount since t=0

    if (TRY_TO_WRITE_NETCDF) {
        const ChemTable* const SM_basis_full = ICS_full_->SM_basis_.get();
        const ChemTable* const SM_all_full = ICS_full_->SM_all_.get();

        // Note: If TRY_TO_WRITE_NETCDF == 1, we only store concentrations of surface complexes and ion exchange species.
        for (int i = 0; i < SM_basis_full->noRows_; ++i) {
            const std::string species_name = SM_basis_full->row_name_[i];
            if (TRY_TO_WRITE_NETCDF >= 2 || (TRY_TO_WRITE_NETCDF == 1 && SM_basis_full->type_[i] != GeochemicalComponentType::AQUEOUS_COMPLEX)) {
                block_concentrations_for_individual_species[species_name] = std::vector<double>(nbl, 0.0);
            }
        }
        for (int i = 0; i < SM_all_full->noRows_; ++i) {
            const std::string species_name = SM_all_full->row_name_[i];
            if (TRY_TO_WRITE_NETCDF >= 2 || (TRY_TO_WRITE_NETCDF == 1 && SM_all_full->type_[i] != GeochemicalComponentType::AQUEOUS_COMPLEX)) {
                block_concentrations_for_individual_species[species_name] = std::vector<double>(nbl, 0.0);
            }
        }

        // We also store precipitated amounts of minerals
        for (int i = 0; i < ICS_full_->size_min_; ++i) {
            const std::string mineral_name = ICS_full_->get_mineral_name(i);
            block_precipitation_for_minerals[mineral_name] = std::vector<double>(nbl, 0.0);
        }

    }

    double pHi = 7.0;
    Ceff.resize(ICS_full_->size_aq_, 0.0);
    scharge.resize(2);
    scharge[0].resize(nbl);
    scharge[1].resize(nbl);

    for (int i = 0; i < ICS_full_->size_aq_; ++i) {
        Cchem_[i] = 0.;
    }
    for (int i = 0; i < ICS_full_->size_aq_; ++i) {
        Ccads_[i] = 0.;
    }

    C_inj.resize(ICS_sol_.size());
    double pH_inj = 0.0;

    const int PRINT_DEBUG = debug;

    for (std::size_t solIdx = 1; solIdx < ICS_sol_.size(); ++solIdx)
    {
        std::string solution_name = ICS_sol_[solIdx]->name();
        std::replace(solution_name.begin(), solution_name.end(), ' ', '_');
        const std::string debug_file = case_name + "_init_" + solution_name;

        // NOTE: We currently enforce charge balance for the initial solutions
        const int chbal_flag = ICS_sol_[solIdx]->charge_balance_;  // backup

        // Make sure that species added by user are available
        C_inj[solIdx].resize(ICS_full_->size_aq_, 0.);
        for (int comp = 0; comp < ICS_sol_[solIdx]->size_aq_; ++comp)
        {
            C_inj[solIdx][ICS_sol_[solIdx]->kcmap_[comp]] = ICS_sol_[solIdx]->c_vchem_[comp];
        }

        std::vector<double> Cads_dum;
        std::vector<double> sch_dum;
        equilibrate_solution(ICS_sol_[solIdx].get(), /* equilibrate=*/ true, debug_file, PRINT_DEBUG, C_inj[solIdx], Cads_dum, pH_inj, sch_dum);

        ICS_sol_[solIdx]->charge_balance_ = chbal_flag;  // set back
    }

    for (int i = 0; i < ICS_full_->size_aq_; ++i) {
        Cchem_[i] = 0.0;
    }
    for (int i = 0; i < ICS_full_->size_aq_; ++i) {
        Ccads_[i] = 0.0;
    }

    for (int i = 0; i < ICS_sol_[0]->size_aq_; ++i)
    {
        assert(ICS_sol_[0]->kcmap_[i] >= 0 && ICS_sol_[0]->kcmap_[i] < static_cast<int>(Cchem_.size()));
        Cchem_[ICS_sol_[0]->kcmap_[i]] = ICS_sol_[0]->c_vchem_[i];
    }

    const double rock_density = ICS_full_->compute_rock_density(porosity);  // kg/m^3

    for (int bufferIdx = 0; bufferIdx < ICS_full_->size_min_; ++bufferIdx)
    {
        const auto name_of_buffer = ICS_full_->get_mineral_name(bufferIdx);
        const double input_concentration = ICS_full_->c_mineral_[bufferIdx];
        const double mol_weight = ICS_full_->SM_mineral_->mol_weight_[bufferIdx];

        if(is_gas_buffer(name_of_buffer))
        {
            Cmin_[bufferIdx] = input_concentration / mol_weight;
        }
        else
        {
            const double mineral_wt_frac = input_concentration;
            const double fac = 1.0e-3*(rock_density/mol_weight)*(1.0-porosity)/porosity;
            Cmin_[bufferIdx] = fac*mineral_wt_frac;
        }
    }

    const double Cw_inv = 1.0;

    // Ion exchange
    const int firstIonExchangerIndex = ICS_full_->size_aq_;
    for (int k = 0; k < ICS_sol_[0]->size_io_; ++k)
    {
        Cchem_[firstIonExchangerIndex + k] = ICS_sol_[0]->c_io_[k] * Cw_inv;
        Ccads_[firstIonExchangerIndex + k] = 0.0;
    }

    // Surface complexes
    const bool include_surface_complexes = ICS_full_->includes_surface_complexes();
    const double specific_surface_area = ICS_full_->SA_;
    const double size_of_diffuse_layer = ICS_full_->d_DL_;

    const int firstSurfaceComplexIndex = ICS_full_->size_aq_ + ICS_full_->size_io_;
    for (int k = 0; k < ICS_sol_[0]->size_surf_; ++k)
    {
        // Convert from sites/nm^2 to mol/L PV
        Cchem_[firstSurfaceComplexIndex + k] = ICS_sol_[0]->c_surf_[k] * (1.0e18 / PhysicalConstants::Avogadro) * specific_surface_area * Cw_inv;
        Ccads_[firstSurfaceComplexIndex + k] = 0.0;
    }

    if (include_surface_complexes)
    {
        sigma_ = scharge[0][0];
        psi_ = 1.0e-3 * scharge[1][0];
    }

    // frac volume of diffusive layer
    double frac_DL = 0.;
    if (size_of_diffuse_layer > 0.0)
    {
        //double Sw = 1.0;
        frac_DL = 1.0e-6 * specific_surface_area * size_of_diffuse_layer / mass_phase_[gFluidPhase::WATER];
    }

    for (int i = 0; i < nbl; ++i)
    {
        Cmin.push_back(Cmin_);  // store minerals in core before calling solver
    }

    double* ptr_Cmin = nullptr;
    if (!Cmin_.empty())
    {
        ptr_Cmin = Cmin_.data();
    }

    SolveChem_I(Cchem_.data(),
                Ccads_.data(),
                ptr_Cmin,
                ICS_full_->log_a_mineral_.data(),
                Temp,
                Pres,
                porosity,
        /* dt= */ 0.0,
                specific_surface_area,
                frac_DL,
                mass_phase_,
                pHi,
                sigma_,
                psi_);

    GCS_->writeCurrentBasVecToFile(case_name + "init.out");
    C_inj[0].resize(ICS_full_->size_aq_, 0.);

    //std::vector<double> mHTO;
    //mHTO.resize(ICS_full_->size_aq_, 0.);
    //mHTO[ICS_full_->get_species_index("H")] = ICS_full_->WHTO_;
    for (int c = 0; c < ICS_full_->size_aq_; ++c) {
        Caq_[c] = Cchem_[c] - Ccads_[c];
        C_inj[0][c] = Caq_[c];

        // At least for now, get directly from solver:
        Cdl_[c] = GCS_->cT_DL_[c];
    }

    // Store solution & mineral amounts
    for (int i = 0; i < nbl; ++i)
    {
        Caq.push_back(Caq_);
        C_DL.push_back(Cdl_);
        Caq_old.push_back(Caq_);
        Ctot.push_back(Cchem_);
        Cads.push_back(Ccads_);
        pH.push_back(pHi);
        scharge[0][i] = sigma_;
        scharge[1][i] = psi_;

        // After calling SolveChem_I(), Cmin_ holds the precipitated amount of minerals...
        for (std::size_t c = 0; c < Cmin_.size(); ++c)
        {
            Cmin[i][c] += Cmin_[c];
        }

        // Store dissolved / precipitated mineral amounts for later output?
        GCS_->Vchem_->update_before_write_solution();
        for (const auto& iter : block_precipitation_for_minerals) {
            const std::string mineral_name = iter.first;
            block_precipitation_for_minerals.at(mineral_name)[i] = GCS_->Vchem_->get_delta_mineral(mineral_name);
        }

    }

    std::vector<double> wksp;
#ifdef NETCDF
    if (TRY_TO_WRITE_NETCDF) {
        assert(!ptrNetCDF_); // file should not already be open...
        ptrNetCDF_ = std::make_unique<NetCDFWriter>(case_name + "_NC");
        std::vector<double> xCoord;
        for (int ix = 0; ix < nbl; ++ix) {
            //Compute midpoints of grid blocks based on input grid block lengths
            //const auto xi = (ix == 0) ? 0.5*dx[0] : x_coord[ix-1] + 0.5*(dx[ix-1] + dx[ix]);
            xCoord.push_back(ix + 1);
        }
        ptrNetCDF_->initialize(nbl, &xCoord[0]);
        wksp = std::vector<double>(nbl, 0.0);
    }
#endif

    std::ofstream off(case_name + "OneDEff.out", std::ios::out);
    std::ofstream off_mass(case_name + "OneDEffMass.out", std::ios::out);
    std::vector<std::string> simulation_species;
    off << "Time\tPVinj\t";
    off_mass << "Time\tPVinj\tVolwat\t";
    for (int i = 0; i < ICS_full_->size_aq_; ++i)
    {
        off << ICS_full_->SM_basis_->row_name_[ICS_full_->pos_[i]] << "\t";
        off_mass << ICS_full_->SM_basis_->row_name_[ICS_full_->pos_[i]] << "\t";
        simulation_species.push_back(ICS_full_->SM_basis_->row_name_[ICS_full_->pos_[i]]);
    }
    off << "pH\n";
    off_mass << "pH\n";

    // Prepare for output

    // A bit hack-ish: We use a lambda to create a function for producing NetCDF output...
    auto write_NetCDF = [this,  // <-- allows changing the state of the current object
        nbl,
        Temp,
        // Local variables which can be modified (i.e., by reference)
        &include_surface_complexes,
        &wksp,
        &pH,
        &Cmin,
        &Caq,
        &Cads,
        &C_DL,
        &block_concentrations_for_individual_species,
        &block_precipitation_for_minerals](double t_hours, double pv_sim)
    {
      if (!ptrNetCDF_) return;

      ptrNetCDF_->incrementTime(t_hours, pv_sim);

      for (int c = 0; c < ICS_full_->size_aq_; ++c) {
          const std::string species_name = ICS_full_->SM_basis_->row_name_[ICS_full_->pos_[c]];
          const std::string species_key = remove_charge_from_species_name(species_name);

          for (auto cellIdx = 0; cellIdx < nbl; ++cellIdx) {
              wksp[cellIdx] = Caq[cellIdx][c];
          }
          ptrNetCDF_->writeField(species_key, &wksp[0]);

          for (auto cellIdx = 0; cellIdx < nbl; ++cellIdx) {
              wksp[cellIdx] = Cads[cellIdx][c];
          }
          ptrNetCDF_->writeField("ads_" + species_key, &wksp[0]);

          if (include_surface_complexes)
          {
              for (auto cellIdx = 0; cellIdx < nbl; ++cellIdx) {
                  wksp[cellIdx] = C_DL[cellIdx][c];
              }
              ptrNetCDF_->writeField("DL_" + species_key, &wksp[0]);
          }
      }

      const double RT_div_F = (PhysicalConstants::IdealGasConstant * Temp) / PhysicalConstants::Faraday;
      for (const auto& iter : block_concentrations_for_individual_species)
      {
          const bool is_potential = (iter.first == "E");
          // Note: If basis species E is present, we print psi too!
          if (is_potential) {
              for (auto cellIdx = 0; cellIdx < nbl; ++cellIdx)
              {
                  const double conc_E0 = iter.second[cellIdx];
                  if(conc_E0 <= 0.0) wksp[cellIdx] = 0.0;
                  else wksp[cellIdx] = 1.0e3 * RT_div_F * std::log(conc_E0);  // mV
              }
              ptrNetCDF_->writeField("psi", &wksp[0]);
          }
          for (auto cellIdx = 0; cellIdx < nbl; ++cellIdx) {
              wksp[cellIdx] = iter.second[cellIdx];
          }
          ptrNetCDF_->writeField("m_" + iter.first, &wksp[0]);
      }

      // We MIGHT wish to print out the net amount of precipitated minerals
      for (const auto& iter : block_precipitation_for_minerals) {
          for (auto cellIdx = 0; cellIdx < nbl; ++cellIdx) {
              wksp[cellIdx] = iter.second[cellIdx];
          }
          ptrNetCDF_->writeField("delta_" + iter.first, &wksp[0]);
      }

      // We always print out the concentration of minerals
      for (auto bufferIndex = 0; bufferIndex < ICS_full_->size_min_; ++bufferIndex) {
          for (auto cellIdx = 0; cellIdx < nbl; ++cellIdx) {
              wksp[cellIdx] = Cmin[cellIdx][bufferIndex];
          }
          ptrNetCDF_->writeField("min_" + to_upper_case(ICS_full_->mineral_name_[bufferIndex]), &wksp[0]);
      }

      ptrNetCDF_->writeField("pH", &pH[0]);
    };

    // Helper function to update concentrations of individual species
    auto update_block_concentrations = [this,
        &block_concentrations_for_individual_species,
        &block_precipitation_for_minerals]
        (std::size_t blockIdx)
    {
      // Before moving on to the next cell, check whether to store concentrations for individual species
      for (const auto& iter : block_concentrations_for_individual_species)
      {
          const std::string species_name = iter.first;
          block_concentrations_for_individual_species.at(species_name)[blockIdx] = GCS_->Vchem_->get_species_concentration(species_name);
      }

      for (const auto& iter : block_precipitation_for_minerals)
      {
          const std::string mineral_name = iter.first;
          // Note: Since we store the cumulative amount of dissolution/precipitation, we ADD to the existing value here
          block_precipitation_for_minerals.at(mineral_name)[blockIdx] += GCS_->Vchem_->get_delta_mineral(mineral_name);
      }
    };

    for(int i=0; i < NoBlocks; ++i) update_block_concentrations(i);
    write_NetCDF(0.0, 0.0);  // t = 0

    // Start simulation loop
    std::vector<EffluentIonData> backup_all_effluent_ion_conc;

    const clock_t begin_time = clock();

    double t_sim = 0.0;
    double pv_sim = 0.0;
    double t_epoch = 0.0;
    int t_step = 0;
    while (!simulation_schedule.simulationHasFinished())
    {
        // Get boundary conditions
        const auto& inlet_boundary_conditions = simulation_schedule.getCurrentEpoch().inlet_bc_;
        const int injected_solution_index = inlet_boundary_conditions.solution_index_;
        assert(injected_solution_index >= 0);
        double* C_in = C_inj[injected_solution_index].data();
        const double current_flow_rate = inlet_boundary_conditions.flow_rate_;

        // Calculate maximum constant time step for advection
        const double residence_time = BlockVolume / current_flow_rate;
        const double dt_advect_max = Flush * residence_time;
        // Possible extension: Allow shut-in periods?

        // Calculate actual time steps (for advection & geochemistry) + increment times
        const double remainining_time_of_epoch = simulation_schedule.getTimeOfCurrentEpoch() - t_epoch;
        assert(remainining_time_of_epoch >= 0);

        double dt_advect = dt_advect_max;
        if (dt_advect >= remainining_time_of_epoch) {
            dt_advect = remainining_time_of_epoch;
            t_epoch = 0.0;
            simulation_schedule.advanceToNextEpoch();
        }
        else {
            t_epoch += dt_advect;
        }

        t_sim += dt_advect;
        pv_sim += current_flow_rate * dt_advect / Volume;
        ++t_step;
        const double dt_solver = Imp ? dt_advect / (1.0 + dt_advect / residence_time) : dt_advect;

        // Do transport + geochemistry, one block at a time!
        for (int i = 0; i < nbl; ++i)
        {
            if (Imp)
            {
                for (int c = 0; c < ICS_full_->size_aq_; ++c)
                {
                    Caq[i][c] = (Caq[i][c] + C_in[c] * Flush) / (1. + Flush);
                    Cads[i][c] /= (1. + Flush);
                }
            }
            else
            {
                for (int c = 0; c < ICS_full_->size_aq_; ++c)
                {
                    Caq_old[i][c] = Caq[i][c];
                    Caq[i][c] = (1. - Flush) * Caq[i][c] + C_in[c] * Flush;
                }
            }

            for (int c = 0; c < ICS_full_->size_io_; ++c)
            {
                if(Imp)
                    Ctot[i][firstIonExchangerIndex + c] = ICS_full_->c_io_[c]/(1+Flush);
                else
                    Ctot[i][firstIonExchangerIndex + c] = ICS_full_->c_io_[c];
            }

            for (int c = 0; c < ICS_full_->size_surf_; ++c)
            {
                if(Imp)
                    Ctot[i][firstSurfaceComplexIndex + c] = ICS_full_->c_surf_[c] * (1.0e18 / PhysicalConstants::Avogadro) * specific_surface_area / (1 + Flush);
                else
                    Ctot[i][firstSurfaceComplexIndex + c] = ICS_full_->c_surf_[c] * (1.0e18 / PhysicalConstants::Avogadro) * specific_surface_area;
            }
            if (include_surface_complexes)
            {
                sigma_ = scharge[0][i];
                psi_ = scharge[1][i] * 1e-3;
            }

            for (int c = 0; c < ICS_full_->size_aq_; ++c) {
                Ctot[i][c] = Caq[i][c] + Cads[i][c];
                if (Ctot[i][c] < NumericalConstants::ZERO_CONCENTRATION_THRESHOLD) {
                    Ctot[i][c] = 0.0;
                }
            }

            std::vector<double> C_dummy;
            if (!Cmin[i].empty())
            {
 //               if (Imp)
 //               {
 //                   for (int c = 0; c < Cmin[i].size(); ++c)
 //                   {
 //                       Cmin[i][c] /= (1. + Flush);
 //                   }
 //               }
                C_dummy = Cmin[i];
                assert(C_dummy.size() == Cmin[i].size());
                ptr_Cmin = &(C_dummy[0]);
            }
            else {
                ptr_Cmin = nullptr;
            }

            SolveChem_I(Ctot[i].data(),
                        Cads[i].data(),
                        ptr_Cmin,
                        ICS_full_->log_a_mineral_.data(),
                        Temp,
                        Pres,
                        porosity,
                        dt_solver,
                        specific_surface_area, //ICS_full_->SA_,
                        frac_DL,
                        mass_phase_,
                        pH[i],
                        sigma_,
                        psi_);

            if (i == punch) // write block
            {
                const std::string file_name = "punch_block_" + std::to_string(i) + "t_step_" + std::to_string(t_step) + ".out";
                GCS_->writeCurrentBasVecToFile(file_name);
            }
            double min_corr = 1.;
            if (Imp) min_corr = (1. + Flush);
            for (std::size_t c = 0; c < C_dummy.size(); ++c)
            {
                Cmin[i][c] += C_dummy[c]*min_corr;
            }

            for (int c = 0; c < ICS_full_->size_aq_; ++c)
            {
                C_DL[i][c] = GCS_->cT_DL_[c]; // At least for now, get directly from solver
                Caq[i][c] = Ctot[i][c] - Cads[i][c];
                if (Imp)
                    Cads[i][c] *= (1 + Flush);
            }

            if (Imp) {
                C_in = &(Caq[i][0]);
            }
            else {
                C_in = &(Caq_old[i][0]);
            }
            update_block_concentrations(i);
        }

        std::cout << "t_sim=" << t_sim / 3600.0 << "hours\n";
        off << t_sim / 3600.0 << "\t";
        off_mass << t_sim / 3600.0 << "\t";
        off << pv_sim << "\t";
        off_mass << pv_sim << "\t";
        double volwat = current_flow_rate * dt_advect;
        off_mass << volwat << "\t";
        for (int c = 0; c < ICS_full_->size_aq_; ++c) {
            off << Caq[nbl - 1][c] << "\t";
            off_mass << Caq[nbl - 1][c] * volwat << "\t";
        }
        off << pH[nbl - 1] << "\n";
        off_mass << pH[nbl - 1] << "\n";

        if (t_sim > 0.0) write_NetCDF(t_sim / 3600.0, pv_sim);

        if(SERIALIZE)
        {
            EffluentIonData tmp_ion_conc;
            for (int c = 0; c < ICS_full_->size_aq_; ++c)
            {
                const auto& ion_name =  ICS_full_->SM_basis_->row_name_[ICS_full_->pos_[c]];
                tmp_ion_conc.names_.push_back(ion_name);
                tmp_ion_conc.values_.push_back(Caq[nbl - 1][c]);
            }
            tmp_ion_conc.names_.push_back("pH");
            tmp_ion_conc.values_.push_back(pH[nbl - 1]);

            backup_all_effluent_ion_conc.push_back(tmp_ion_conc);
        }
    } // end simulation loop

    std::cout << "Time used " << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";
    std::cout << "Number of calls to solver: " << GCS_->noSolverCalls() << "(" << nbl << " blocks)\n";

    solutions_manager_.print_information();

    if (SERIALIZE >= 2)
    {
        // Store contents of BasVecInfo in a special file to be able to reconstruct it later (for testing)
        std::ofstream cereal_file(case_name + "_transport.cereal", std::ios::binary);
        cereal::BinaryOutputArchive archive(cereal_file);
        archive(backup_all_effluent_ion_conc);
    }
    return backup_all_effluent_ion_conc;
}

std::string ThermoTableSolver::solve(const std::string& case_name,
                                     std::istream& inputStream)
{
    (void)case_name;

    std::map<std::string, std::string> simple_key_value_pairs;
    simple_key_value_pairs["Species"] = "H2O";
    simple_key_value_pairs["Temp"] = "25.0";
    simple_key_value_pairs["Pres"] = "1.0e5";
    simple_key_value_pairs["TempUnit"] = "C";
    simple_key_value_pairs["PresUnit"] = "Pa";
    simple_key_value_pairs["IncludeIndex"] = "1";

    InputReader input_reader(simple_key_value_pairs, false);
    input_reader.define_block_keyword("SPECIESLIST", {"/ END", "/END"});
    input_reader.define_block_keyword("TEMPS", {"/ END", "/END"});
    input_reader.define_block_keyword("PRESSURES", {"/ END", "/END"});
    input_reader.read(inputStream);

    auto species_names = readSpeciesNames(input_reader);
    auto temperatures_celsius = readScalarOrBlock(input_reader, "Temp", "TEMPS");
    auto pressures_bar = readScalarOrBlock(input_reader, "Pres", "PRESSURES");

    convertTemperaturesToCelsius(temperatures_celsius, input_reader.get_simple_keyword_value("TempUnit"));
    convertPressuresToBar(pressures_bar, input_reader.get_simple_keyword_value("PresUnit"));

    const bool include_index = parseBoolKeyword(input_reader.get_simple_keyword_value("IncludeIndex"), "IncludeIndex");

    ThermoTableCalculator calculator;
    std::ostringstream output;
    for (std::size_t i = 0; i < species_names.size(); ++i)
    {
        const auto rows = calculator.evaluate(species_names[i], temperatures_celsius, pressures_bar);

        if (i > 0)
        {
            output << '\n';
        }

        output << "SPECIES " << species_names[i] << '\n';
        output << ThermoTableCalculator::formatRows(rows, include_index);

        if (ThermoTableCalculator::isWaterSpeciesName(species_names[i]))
        {
            std::vector<double> row_temperatures;
            row_temperatures.reserve(rows.size());
            for (const auto& row : rows)
            {
                row_temperatures.push_back(row.T);
            }

            output << "IAPWS97_PSAT_H2O\n";
            output << calculator.evaluateWaterSaturationPressureFormatted(row_temperatures, include_index);
        }
    }

    return output.str();
}

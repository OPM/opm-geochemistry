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
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/geochemistry/OpmGeoChemInterface.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>
#include <opm/simulators/geochemistry/Core/GeoChemPhases.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <stdexcept>

namespace {

// Print warning helper
void printWarning (const std::string& file_name, const std::string& block, const char* keyword, bool ignore)
{
    std::string msg;
    if (ignore) {
        msg = fmt::format("Block \"{}\" in JSON file {} will be ignored! Use keyword {} instead!",
                            block, file_name, keyword);
    }
    else {
        msg = fmt::format("Block \"{}\" in JSON file {} will be overwritten by keyword {}!",
                            block, file_name, keyword);
    }
    Opm::OpmLog::warning(msg);
}

}

// ////
// PUBLIC METHODS
// ///
void OpmGeoChemInterface::initialize_from_opm_deck(const std::string& file_name,
                                                   const std::vector<std::string>& species,
                                                   const std::vector<std::string>& minerals,
                                                   const std::vector<std::string>& ion_ex,
                                                   bool charge_balance,
                                                   std::pair<double, double> tol,
                                                   int splay_tree_resolution)
{
    // Append species from OPM deck to JSON
    nlohmann::json appended_json = opmDeckSpeciesToJSON_(file_name, charge_balance, species, minerals, ion_ex);

    // Parse JSON
    const std::string json_string = appended_json.dump();
    auto json_parsed = geoChemPhases_.resetFromJson(json_string);
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const int max_size = geoChemPhases_.selectPhases(geoChemPhases_.getPhaseNamesAndTypes(), phases_to_be_used);
    set_geochemical_database_modifications_from_json(json_string);

    // Init.
    ICS_full_ = InitChem::create_from_input_data(nullptr,
                                                 db_changes_made_by_user_,
                                                 phases_to_be_used,
                                                 max_size,
                                                 "GEOCHEM",
                                                 species,
                                                 splay_tree_resolution);
    allocate_memory_for_solver_and_splay_tree_etc();

    // Solver tolerances: tol = {mbal, ph}
    if (json_parsed.contains("chemtol") || json_parsed.contains("CHEMTOL")) {
        printWarning(file_name, "CHEMTOL", "GEOCHEM", /*ignore=*/false);
    }
    set_solver_tolerances(*GCS_, tol.first, tol.second);

    // Check species defined in OPM deck vs. geochemical solver
    checkUserOrder_(species);
}

void OpmGeoChemInterface::initialize_json(const std::string& file_name,
                                          const std::vector<std::string>& user_order)
{
    // Read input from JSON file
    auto json_parsed = geoChemPhases_.resetFromJson(file_name);
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const std::map<std::string, std::vector<std::string>> db_changes_made_by_user;
    const int max_size = geoChemPhases_.selectPhases(geoChemPhases_.getPhaseNamesAndTypes(), phases_to_be_used);

    // Init
    ICS_full_ = InitChem::create_from_input_data(nullptr,
                                                 db_changes_made_by_user,
                                                 phases_to_be_used,
                                                 max_size,
                                                 "GEOCHEM",
                                                 user_order,
                                                 0);
    allocate_memory_for_solver_and_splay_tree_etc();

    // Solver option
    // TODO: Add to JSON parser instead of function below?
    // Add option to change ICS_full_->INTERPOLATE ?
    modify_selected_solver_options_json(*GCS_, json_parsed);  // HACK(?)

    // Check species defined in OPM deck vs. geochemical solver
    checkUserOrder_(user_order, file_name);
}

void OpmGeoChemInterface::initialize(const std::string& file_name,
                                     double temperature,
                                     double porosity,
                                     const std::vector<std::string>& user_order)
{
    if (ICS_full_) throw MultipleInitializationException("Cannot initialize ICS_full_, it already exists...");

    std::ifstream inputStream(file_name, std::ios::binary);
    ptrInputReader_->read(inputStream);
    set_geochemical_database_modifications_from_user_input();

    geoChemPhases_.resetFromInputStream(inputStream);

    // Get all solutions entered, and mix all solutions into one solution (same for minerals and surfaces)
    // NB: Must set interpolate flag before allocating memory for solver, splay tree, etc.
    create_ICS_full(user_order);

    ICS_full_->INTERPOLATE_ = std::stoi(ptrInputReader_->get_simple_keyword_value("INTERPOLATE"));
    ICS_full_->PRINT_DEBUG_CHEM_ = std::stoi(ptrInputReader_->get_simple_keyword_value("DEBUG"));
    ICS_full_->possibly_change_inconsistent_options();

    allocate_memory_for_solver_and_splay_tree_etc();
    modify_selected_solver_options(*GCS_);

    Temp_ = temperature; // ICS_full_->Temp_;

    // We start by adding all phases initially in the reservoir, then we add the injected solutions.
    auto [reservoir_phases, names_of_injected_solutions] = geoChemPhases_.GetInjectReservoirSimple();
    AddPhases(reservoir_phases, "RESERVOIR");
    for (const auto& solution_name: names_of_injected_solutions)
    {
        AddAqueousSolution(solution_name);
    }

    for (int solutionIndex = 0; solutionIndex < static_cast<int>(ICS_sol_.size()); ++solutionIndex)
    {
        // Equilibration is done locally in each block (later).
        set_data_for_solution(solutionIndex, porosity);
    }
}

std::vector<std::size_t>
OpmGeoChemInterface::mineral_index_map(const std::vector<std::string>& deck_minerals) const
{
    std::vector<std::size_t> map;
    map.reserve(deck_minerals.size());
    std::vector<bool> used(ICS_full_->size_min_, false);
    for (const auto& name : deck_minerals) {
        const int idx = ICS_full_->get_mineral_index(name);
        if (idx < 0) {
            const std::string msg
                = fmt::format("MINERAL: {} not found in the geochemistry solver.", name);
            Opm::OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        if (used[idx]) {
            const std::string msg = fmt::format("MINERAL: duplicate {}.", name);
            Opm::OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        used[idx] = true;
        map.push_back(static_cast<std::size_t>(idx));
    }
    return map;
}

void OpmGeoChemInterface::calculate_initial_mineral_concentration(std::vector<double>& Cmin,
                                                                  double porosity,
                                                                  const std::unordered_map<std::string, double>& weight_mineral)
{
    if (weight_mineral.size() != static_cast<std::size_t>(ICS_full_->size_min_)) {
        const std::string msg
            = fmt::format("Weight fractions are given for {} minerals, but the geochemistry solver "
                          "has {}.",
                          weight_mineral.size(),
                          ICS_full_->size_min_);
        Opm::OpmLog::error(msg);
        throw std::runtime_error(msg);
    }

    for (const auto& [name, val] : weight_mineral) {
        int minIdx = ICS_full_->get_mineral_index(name);
        if (minIdx < 0) {
            const std::string msg = fmt::format("Mineral = {} not found internally in geochemistry solver!", name);
            Opm::OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        ICS_full_->c_mineral_[minIdx] = val;
    }

    const double rock_density = ICS_full_->compute_rock_density(porosity);  // kg/m^3

    Cmin.resize(ICS_full_->size_min_);
    for (int bufferIdx = 0; bufferIdx < ICS_full_->size_min_; ++bufferIdx)
    {
        const auto name_of_buffer = ICS_full_->get_mineral_name(bufferIdx);
        const double input_concentration = ICS_full_->c_mineral_[bufferIdx];
        const double mol_weight = ICS_full_->SM_mineral_->mol_weight_[ICS_full_->pos_min_[bufferIdx]];

        if(is_gas_buffer(name_of_buffer))
        {
            Cmin[bufferIdx] = input_concentration / mol_weight;
        }
        else
        {
            const double mineral_wt_frac = input_concentration;
            const double fac = 1.0e-3*(rock_density/mol_weight)*(1.0-porosity)/porosity;
            Cmin[bufferIdx] = fac*mineral_wt_frac;
        }
    }
}

void OpmGeoChemInterface::set_surface_concentrations(double swat,
                                                     std::vector<double>& C_tot,
                                                     double& frac_DL,
                                                     const std::unordered_map<std::string, double>& C_io)
{
    // Set ion exchange concentration
    if (C_io.size() != static_cast<std::size_t>(ICS_full_->size_io_)) {
        const std::string msg
            = fmt::format("Concentrations are given for {} ion exchange species, but the "
                          "geochemistry solver has {}.",
                          C_io.size(),
                          ICS_full_->size_io_);
        Opm::OpmLog::error(msg);
        throw std::runtime_error(msg);
    }

    for (int i = 0; i < ICS_full_->size_io_; ++i) {
        const auto& name = ICS_full_->io_name_[i];
        ICS_full_->c_io_[i] = C_io.at(name);
    }

    return SetSurfaceConc(swat, C_tot, frac_DL);
}

void OpmGeoChemInterface::set_solver_tolerances(GCSolver& GCS_in,
                                                double mbal_tol,
                                                double ph_tol)
{
    GCS_in.options_.MBAL_CONV_CRITERION_ = mbal_tol;
    GCS_in.options_.PH_CONV_CRITERION_ = ph_tol;
}

std::vector<double>& OpmGeoChemInterface::get_log_a_mineral()
{ return ICS_full_->log_a_mineral_; }

const std::vector<double>& OpmGeoChemInterface::get_log_a_mineral() const
{ return ICS_full_->log_a_mineral_; }

// ////
// PRIVATE METHODS
// ///
nlohmann::json OpmGeoChemInterface::opmDeckSpeciesToJSON_(const std::string& file_name,
                                                          bool charge_balance,
                                                          const std::vector<std::string>& species,
                                                          const std::vector<std::string>& minerals,
                                                          const std::vector<std::string>& ion_ex)
{
    // SPECIES keyword required!
    if (species.empty()) {
        const std::string msg = "SPECIES keyword required in geochemistry solver!";
        Opm::OpmLog::error(msg);
        throw std::runtime_error(msg);
    }

    nlohmann::json json_opm;
    if (!file_name.empty()) {
        // Read JSON file
        std::ifstream in(file_name);
        if (!in) {
            const std::string msg = fmt::format("Could not open JSON file: {}", file_name);
            Opm::OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        in >> json_opm;
    }

    // Append SOLUTION 0 species (i.e. aqueous species)
    const std::string solution_block = PhaseKeyword(GeochemicalPhaseType::AQUEOUS_SOLUTION);
    if (json_opm.contains(solution_block) || json_opm.contains(to_lower_case(solution_block))) {
        printWarning(file_name, solution_block, "SPECIES", /*ignore=*/false);
    }
    json_opm[solution_block] = nlohmann::json::object();
    json_opm[solution_block][solution_block + " 0"] = nlohmann::json::object();
    for (const auto& elem : species) {
        if (elem == "H" && charge_balance) {
            json_opm[solution_block][solution_block + " 0"][elem] = "1.0 charge";
        }
        else {
            json_opm[solution_block][solution_block + " 0"][elem] = "1e-12";  // arbitrary value!
        }
    }

    // Append minerals
    const std::string mineral_block = PhaseKeyword(GeochemicalPhaseType::EQUILIBRIUM_MINERAL);
    if (!minerals.empty()) {
        if (json_opm.contains(mineral_block) || json_opm.contains(to_lower_case(mineral_block))) {
            printWarning(file_name, mineral_block, "MINERAL", /*ignore=*/false);
        }
        json_opm[mineral_block] = nlohmann::json::object();
        json_opm[mineral_block][mineral_block + " 0"] = nlohmann::json::object();
        for (const auto& elem : minerals) {
            json_opm[mineral_block][mineral_block + " 0"][elem] = "1";  // arbitrary value!
        }
    }
    else {
        if (json_opm.contains(mineral_block) || json_opm.contains(to_lower_case(mineral_block))) {
            printWarning(file_name, mineral_block, "MINERAL", /*ignore=*/true);
            json_opm.erase(mineral_block);
        }
    }

    // Append ion exchange
    const std::string iexchange_block = PhaseKeyword(GeochemicalPhaseType::EXCHANGE_SITES);
    if (!ion_ex.empty()) {
        if (json_opm.contains(iexchange_block) || json_opm.contains(to_lower_case(iexchange_block))) {
            printWarning(file_name, iexchange_block, "IONEX", /*ignore=*/false);
        }
        json_opm[iexchange_block] = nlohmann::json::object();
        json_opm[iexchange_block][iexchange_block + " 0"] = nlohmann::json::object();
        for (const auto& elem : ion_ex) {
            json_opm[iexchange_block][iexchange_block + " 0"][elem] = "";  // sets default value
        }
    }
    else {
        // Delete ion exchange from JSON
        if (json_opm.contains(iexchange_block)|| json_opm.contains(to_lower_case(iexchange_block))) {
            printWarning(file_name, iexchange_block, "IONEX", /*ignore=*/true);
            json_opm.erase(iexchange_block);
        }
    }

    return json_opm;
}

void OpmGeoChemInterface::checkUserOrder_(const std::vector<std::string>& user_order,
                                          const std::optional<std::string>& file_name)
{
    const auto& basis = *ICS_full_->SM_basis_;

    auto join = [](const std::vector<std::string>& names) {
        std::string joined;
        for (const auto& name : names) {
            joined += (joined.empty() ? "" : " ") + name;
        }
        return joined;
    };

    // Names used in the messages. In the OPM deck the species are given by the SPECIES keyword and
    // the geochemistry solver is built from the deck. When a JSON input is given, the species list
    // is passed in separately and the geochemistry solver's species come from the JSON input. This
    // is either the path of a file, which is shown, or the JSON text itself, which can be long.
    const bool has_file = file_name.has_value();
    const std::string species_list = has_file ? "Species list" : "SPECIES";
    const std::string solver_source = [&]() -> std::string {
        if (!has_file) {
            return "the geochemistry solver";
        }
        const auto first = file_name->find_first_not_of(" \t\r\n");
        const bool is_json_text = first != std::string::npos && (*file_name)[first] == '{';
        return is_json_text ? "the JSON input" : fmt::format("\"{}\"", *file_name);
    }();

    // Database rows of the species in SPECIES. Names are looked up in the database, so a nickname
    // and the full name of the same species (e.g. NA and Na+), or different capitalisation, give
    // the same row. A species given more than once is an error in SPECIES that the geochemistry
    // solver cannot resolve: it either ends up with a different number of species or with one
    // species twice.
    std::vector<std::string> duplicates;
    std::vector<int> user_rows;
    user_rows.reserve(user_order.size());
    for (const auto& name : user_order) {
        const int row = basis.get_row_index(name);
        const auto earlier = std::ranges::find(user_rows, row);
        if (row >= 0 && earlier != user_rows.end()) {
            duplicates.push_back(fmt::format(
                "{} (same species as {})", name, user_order[earlier - user_rows.begin()]));
        }
        user_rows.push_back(row);
    }

    // Database rows of the species of the geochemistry solver, and the species that it has but
    // SPECIES does not (e.g. required by a mineral)
    std::vector<int> solver_rows;
    solver_rows.reserve(ICS_full_->size_aq_);
    std::vector<std::string> added_by_solver;
    for (int i = 0; i < ICS_full_->size_aq_; ++i) {
        const int solver_row = basis.get_row_index(ICS_full_->basis_species_name_[i]);
        solver_rows.push_back(solver_row);
        if (std::ranges::find(user_rows, solver_row) == user_rows.end()) {
            added_by_solver.push_back(ICS_full_->basis_species_name_[i]);
        }
    }

    // Species in SPECIES that the geochemistry solver does not have: either unknown in the
    // database, or known but not among the species of the geochemistry solver
    std::vector<std::string> not_in_solver;
    for (std::size_t i = 0; i < user_order.size(); ++i) {
        if (user_rows[i] < 0 || std::ranges::find(solver_rows, user_rows[i]) == solver_rows.end()) {
            not_in_solver.push_back(user_order[i]);
        }
    }

    // The species must be the same as in the geochemistry solver, and in the same order as the
    // concentrations are mapped to the geochemistry solver species index. Everything that is
    // wrong is reported at once.
    const bool same_species = duplicates.empty() && added_by_solver.empty() && not_in_solver.empty();
    if (same_species && user_rows == solver_rows) {
        return;
    }

    std::vector<std::string> lines;
    if (!duplicates.empty()) {
        lines.push_back(fmt::format("{}: duplicate {}.", species_list, join(duplicates)));
    }
    if (!not_in_solver.empty()) {
        lines.push_back(fmt::format(
            "{}: remove {}, not in {}.", species_list, join(not_in_solver), solver_source));
    }
    if (!added_by_solver.empty()) {
        lines.push_back(fmt::format(
            "{}: add {}, required by {}.", species_list, join(added_by_solver), solver_source));
    }
    if (same_species) {
        const std::vector<std::string> solver_species(ICS_full_->basis_species_name_.begin(),
                                                      ICS_full_->basis_species_name_.begin()
                                                          + ICS_full_->size_aq_);
        lines.push_back(fmt::format("{}: order differs from {}: {} vs {}.",
                                    species_list,
                                    solver_source,
                                    join(user_order),
                                    join(solver_species)));
    }
    std::string msg;
    for (const auto& line : lines) {
        msg += (msg.empty() ? "" : "\n") + line;
    }
    Opm::OpmLog::error(msg);
    throw std::runtime_error(msg);
}

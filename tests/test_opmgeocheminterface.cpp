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
#include "config.h"

#define BOOST_TEST_MODULE OpmGeochemInterfaceTests

#include <boost/test/unit_test.hpp>

#include <opm/simulators/geochemistry/OpmGeoChemInterface.hpp>

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


namespace
{

// Initialize an interface from deck-style species lists (no input file, no simulator needed)
std::shared_ptr<OpmGeoChemInterface>
make_interface(const std::vector<std::string>& species,
               const std::vector<std::string>& minerals = {},
               const std::vector<std::string>& ion_ex = {})
{
    auto interface = std::make_shared<OpmGeoChemInterface>();
    interface->initialize_from_opm_deck(/*file_name=*/"",
                                        species,
                                        minerals,
                                        ion_ex,
                                        /*charge_balance=*/true,
                                        std::make_pair(1e-5, 1e-6),
                                        /*splay_tree_resolution=*/10);
    return interface;
}

// Predicate for BOOST_CHECK_EXCEPTION: the exception message must contain the given text
auto
message_contains(const std::string& text)
{
    return [text](const std::exception& e) {
        return std::string(e.what()).find(text) != std::string::npos;
    };
}

// Same, but all the given texts must be in the message
auto
message_contains_all(const std::vector<std::string>& texts)
{
    return [texts](const std::exception& e) {
        const std::string message(e.what());
        return std::ranges::all_of(texts, [&message](const std::string& text) {
            return message.find(text) != std::string::npos;
        });
    };
}

const std::vector<std::string> default_species = {"H", "NA", "SO4", "K", "CL", "CA", "HCO3"};

// JSON input that defines the species H, NA, K and CL for the solver
const std::string json_species = R"({"SOLUTION": {"SOLUTION 0": )"
                                 R"({"pH": "7 charge", "Na": "1e-3", "K": "1e-3", "Cl": "2e-3"}}})";

} // anonymous namespace


BOOST_AUTO_TEST_CASE(InitializeFromDeckTest)
{
    // Setup
    std::string empty_file_name = "";
    bool charge_balance = true;
    const std::pair<double, double> tol = std::make_pair<double, double>(1e-5, 1e-6);
    int splay_tree_resolution = 10;

    // Initialize interface
    std::shared_ptr<OpmGeoChemInterface> geoChemInterface = std::make_shared<OpmGeoChemInterface>();
    const std::vector<std::string> species = { "H", "NA", "SO4", "K", "CL", "CA", "HCO3" };
    const std::vector<std::string> minerals = { "CAL" };
    const std::vector<std::string> ion_ex = { "X" };
    geoChemInterface->initialize_from_opm_deck(empty_file_name,
                                               species,
                                               minerals,
                                               ion_ex,
                                               charge_balance,
                                               tol,
                                               splay_tree_resolution);

    // Check number of species, minerals, and ion exchangers
    BOOST_CHECK_EQUAL(species.size(), geoChemInterface->numberOfAqueousBasisSpecies());
    BOOST_CHECK_EQUAL(minerals.size(), geoChemInterface->numberOfMinerals());
    BOOST_CHECK_EQUAL(ion_ex.size(), geoChemInterface->numberOfIonExchange());
    BOOST_CHECK_EQUAL(species.size() + ion_ex.size(), geoChemInterface->numberOfBasisSpecies());

    // Check if species names are located in interface
    const auto& basis_species = geoChemInterface->GetBasisSpeciesNames();
    for (const auto& elem : species) {
        BOOST_CHECK_MESSAGE(
            std::find(basis_species.begin(), basis_species.end(), elem) != basis_species.end(),
            "Species missing in GetBasisSpeciesNames() = " << elem
        );
    }

    // Check mineral name
    const auto& mineral_species = geoChemInterface->GetMineralNames();
    BOOST_CHECK_MESSAGE(
        std::find(mineral_species.begin(), mineral_species.end(), minerals[0]) != mineral_species.end(),
        "Mineral missing in GetMineralNames() = " << minerals[0]
    );

    // Check ion exchange name
    const auto& ionex_species = geoChemInterface->GetIONames();
    BOOST_CHECK_MESSAGE(
        std::find(basis_species.begin(), basis_species.end(), ion_ex[0]) != basis_species.end(),
        "Ion exchange missing in GetBasisSpeciesNames() = " << ion_ex[0]
    );
    BOOST_CHECK_MESSAGE(
        std::find(ionex_species.begin(), ionex_species.end(), ion_ex[0]) != ionex_species.end(),
        "Ion exchange missing in GetIONames() = " << ion_ex[0]
    );

    // Calculate initial mineral concentration test
    // Formula = total_rock_density / mole_weight_mineral * (1-poro) / poro * weight_fraction_mineral
    std::unordered_map<std::string, double> weight_mineral{ { minerals[0], 0.001 } };
    std::vector<double> Cmin;
    double porosity = 0.25;
    geoChemInterface->calculate_initial_mineral_concentration(Cmin, porosity, weight_mineral);
    double Cmin_expected = 0.08093;
    BOOST_CHECK_CLOSE(Cmin_expected, Cmin[0], 1e-2);

    // Calculate surface concentrations
    // NOTE: ion exchangers are inserted in Ctot _after_ aqueous species
    std::unordered_map<std::string, double> Cion { { ion_ex[0], 1.1e-3 } };
    double swat = 0.5;
    std::size_t nBasis = geoChemInterface->numberOfBasisSpecies();
    std::vector<double> Ctot(nBasis, 0.0);
    double frac_dl = 0.0;
    geoChemInterface->set_surface_concentrations(swat, Ctot, frac_dl, Cion);
    std::size_t nAq = nBasis - 1;
    BOOST_CHECK_CLOSE(Ctot[nAq], (1.0 / swat) * Cion.at(ion_ex[0]), 1e-6);

    // Check get functions
    for (const auto& elem : geoChemInterface->get_log_a_mineral()) {
        BOOST_CHECK_EQUAL(elem, 0.0);
    }
}


BOOST_AUTO_TEST_CASE(InitializeWithoutMineralsAndIonExchangeTest)
{
    const auto interface = make_interface(default_species);

    BOOST_CHECK_EQUAL(interface->numberOfAqueousBasisSpecies(), default_species.size());
    BOOST_CHECK_EQUAL(interface->numberOfMinerals(), 0);
    BOOST_CHECK_EQUAL(interface->numberOfIonExchange(), 0);
    BOOST_CHECK_EQUAL(interface->numberOfBasisSpecies(), default_species.size());

    // The log(a) vector is allocated with the solver's buffer capacity, which can exceed the number
    // of minerals
    BOOST_CHECK_GE(interface->get_log_a_mineral().size(), interface->numberOfMinerals());
    for (const auto elem : interface->get_log_a_mineral()) {
        BOOST_CHECK_EQUAL(elem, 0.0);
    }
}


BOOST_AUTO_TEST_CASE(SpeciesOrderFollowsUserOrderTest)
{
    // The solver must order aqueous species the same way as the OPM deck
    const std::vector<std::string> reordered = {"HCO3", "CA", "H", "CL", "K", "SO4", "NA"};
    const auto interface = make_interface(reordered);

    const auto& basis_species = interface->GetBasisSpeciesNames();
    BOOST_REQUIRE_GE(basis_species.size(), reordered.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(basis_species.begin(),
                                  basis_species.begin() + reordered.size(),
                                  reordered.begin(),
                                  reordered.end());
}


BOOST_AUTO_TEST_CASE(SpeciesOrderAcceptsFullNamesAndNicknamesTest)
{
    // A species can be given by its database name or by its nickname, in any order. The solver
    // keeps the user's order and spelling (in upper case).
    const std::vector<std::string> mixed = {"HCO3", "Ca+2", "H", "Cl-", "K", "SO4", "Na+"};
    const std::vector<std::string> expected = {"HCO3", "CA+2", "H", "CL-", "K", "SO4", "NA+"};
    const auto interface = make_interface(mixed);

    const auto& basis_species = interface->GetBasisSpeciesNames();
    BOOST_REQUIRE_GE(basis_species.size(), expected.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(basis_species.begin(),
                                  basis_species.begin() + expected.size(),
                                  expected.begin(),
                                  expected.end());
}


BOOST_AUTO_TEST_CASE(DuplicateSpeciesThrowsTest)
{
    // The same name twice must be reported as such
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_from_opm_deck("",
                                           {"H", "NA", "NA", "K", "CL", "CA", "HCO3"},
                                           {},
                                           {},
                                           true,
                                           std::make_pair(1e-5, 1e-6),
                                           10),
        std::runtime_error,
        message_contains_all({"duplicate", "NA"}));

    // The same with JSON input, where the species list is given separately
    OpmGeoChemInterface json_interface;
    BOOST_CHECK_EXCEPTION(
        json_interface.initialize_json(json_species, {"H", "NA", "NA", "K", "CL"}),
        std::runtime_error,
        message_contains_all({"Species list", "duplicate", "NA"}));
}


BOOST_AUTO_TEST_CASE(SameSpeciesWithDifferentCapitalisationThrowsTest)
{
    // Na and NA are the same species. The deck parser compares names case sensitively, so this
    // reaches the solver, where the two names collapse into one species.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_from_opm_deck(
            "", {"H", "NA", "Na", "K", "CL", "CA"}, {}, {}, true, std::make_pair(1e-5, 1e-6), 10),
        std::runtime_error,
        message_contains_all({"duplicate", "Na"}));

    // The same with JSON input, where the species list is given separately
    OpmGeoChemInterface json_interface;
    BOOST_CHECK_EXCEPTION(
        json_interface.initialize_json(json_species, {"H", "NA", "Na", "K", "CL"}),
        std::runtime_error,
        message_contains_all({"Species list", "duplicate", "Na"}));
}


BOOST_AUTO_TEST_CASE(SameSpeciesWithNicknameAndFullNameThrowsTest)
{
    // NA is the nickname of Na+. Both are accepted by the deck parser, and without this check the
    // solver ends up with the same database species twice.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_from_opm_deck(
            "", {"H", "NA", "NA+", "K", "CL", "CA"}, {}, {}, true, std::make_pair(1e-5, 1e-6), 10),
        std::runtime_error,
        message_contains_all({"duplicate", "NA+"}));

    // The same with JSON input, where the species list is given separately
    OpmGeoChemInterface json_interface;
    BOOST_CHECK_EXCEPTION(
        json_interface.initialize_json(json_species, {"H", "NA", "NA+", "K", "CL"}),
        std::runtime_error,
        message_contains_all({"Species list", "duplicate", "NA+"}));
}


BOOST_AUTO_TEST_CASE(DuplicateSpeciesWithSolverAddedSpeciesThrowsTest)
{
    // The duplicate NA collapses to one species and the solver adds HCO3 for the mineral CAL, so
    // the number of species equals the number of names given. The solver then keeps its own order
    // (CA CL H K NA HCO3-) instead of the order in SPECIES, and this must be detected.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_from_opm_deck("",
                                           {"H", "NA", "NA", "K", "CL", "CA"},
                                           {"CAL"},
                                           {},
                                           true,
                                           std::make_pair(1e-5, 1e-6),
                                           10),
        std::runtime_error,
        message_contains_all({"SPECIES: duplicate", "NA", "SPECIES: add", "HCO3"}));

    // The same with JSON input, where CL is the species that the JSON input has and the list lacks
    OpmGeoChemInterface json_interface;
    BOOST_CHECK_EXCEPTION(json_interface.initialize_json(json_species, {"H", "NA", "NA", "K"}),
                          std::runtime_error,
                          message_contains_all({"Species list: duplicate",
                                                "NA",
                                                "Species list: add",
                                                "required by the JSON input",
                                                "CL"}));
}


BOOST_AUTO_TEST_CASE(MissingSpeciesRequiredByMineralThrowsTest)
{
    // Calcite (CAL) needs HCO3, which the solver adds itself when it is missing from SPECIES. The
    // error must tell which species to add.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_from_opm_deck(
            "", {"H", "NA", "K", "CL", "CA"}, {"CAL"}, {}, true, std::make_pair(1e-5, 1e-6), 10),
        std::runtime_error,
        message_contains_all({"SPECIES: add", "HCO3"}));

    // The same with JSON input: the list lacks the species K and CL that the JSON input defines
    OpmGeoChemInterface json_interface;
    BOOST_CHECK_EXCEPTION(json_interface.initialize_json(json_species, {"H", "NA"}),
                          std::runtime_error,
                          message_contains_all({"Species list: add", "K CL"}));
}


BOOST_AUTO_TEST_CASE(SpeciesNotInSolverSpeciesThrowsTest)
{
    // CA is a valid species in the database, but it is not among the species that the JSON input
    // defines (H, NA, K, CL). The error must name it as a species to remove.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(interface.initialize_json(json_species, {"H", "NA", "K", "CL", "CA"}),
                          std::runtime_error,
                          message_contains_all({"Species list: remove", "CA"}));
}


BOOST_AUTO_TEST_CASE(MissingSpeciesGivenByFullNameTest)
{
    // Na+ is the full name of NA. Only CL is missing from the species list, so NA must not be
    // reported as missing.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_json(json_species, {"H", "Na+", "K"}),
        std::runtime_error,
        [](const std::exception& e) {
            const std::string message(e.what());
            const auto added = message.find("Species list: add ");
            return added != std::string::npos
                && message.find("CL", added) != std::string::npos
                && message.find("NA", added) == std::string::npos;
        });
}


BOOST_AUTO_TEST_CASE(AllSpeciesErrorsAreReportedTogetherTest)
{
    // The JSON input defines H, NA, K and CL. The list has NA twice, CA that the JSON input does
    // not have, and lacks CL. All three problems must be in the same error.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_json(json_species, {"H", "NA", "NA", "K", "CA"}),
        std::runtime_error,
        [](const std::exception& e) {
            const std::string message(e.what());
            const auto removed = message.find("Species list: remove ");
            const auto added = message.find("Species list: add ");
            return message.find("Species list: duplicate") != std::string::npos
                && removed != std::string::npos
                && message.find("CA", removed) != std::string::npos
                && added != std::string::npos
                && message.find("CL", added) != std::string::npos;
        });
}


BOOST_AUTO_TEST_CASE(SpeciesOrderDiffersThrowsTest)
{
    // The same species as the JSON input (H, NA, K, CL), given by full names in another order than
    // the solver has. Only nicknames are reordered by the solver, so the order is wrong and must
    // be reported, with both lists, and without any species to add or remove.
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_json(json_species, {"Na+", "K+", "H+", "Cl-"}),
        std::runtime_error,
        [](const std::exception& e) {
            const std::string message(e.what());
            return message.find("Species list: order differs from the JSON input: Na+ K+ H+ Cl- vs ")
                       != std::string::npos
                && message.find("add") == std::string::npos
                && message.find("remove") == std::string::npos;
        });
}


BOOST_AUTO_TEST_CASE(InitializeWithoutSpeciesThrowsTest)
{
    // The SPECIES keyword is required
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(
        interface.initialize_from_opm_deck("", {}, {}, {}, true, std::make_pair(1e-5, 1e-6), 10),
        std::runtime_error,
        message_contains("SPECIES keyword required"));
}


BOOST_AUTO_TEST_CASE(InitializeWithMissingJsonFileThrowsTest)
{
    OpmGeoChemInterface interface;
    BOOST_CHECK_EXCEPTION(interface.initialize_from_opm_deck("this_file_does_not_exist.json",
                                                             default_species,
                                                             {},
                                                             {},
                                                             true,
                                                             std::make_pair(1e-5, 1e-6),
                                                             10),
                          std::runtime_error,
                          message_contains("Could not open JSON file"));
}


BOOST_AUTO_TEST_CASE(LogAMineralConstOverloadTest)
{
    const auto interface = make_interface(default_species, {"CAL"});

    // The const overload refers to the same storage as the mutable one
    const OpmGeoChemInterface& const_interface = *interface;
    BOOST_CHECK_EQUAL(&const_interface.get_log_a_mineral(), &interface->get_log_a_mineral());
    BOOST_CHECK_GE(const_interface.get_log_a_mineral().size(), interface->numberOfMinerals());

    // Writes through the mutable overload are visible through the const one
    interface->get_log_a_mineral()[0] = -1.5;
    BOOST_CHECK_EQUAL(const_interface.get_log_a_mineral()[0], -1.5);
}


BOOST_AUTO_TEST_CASE(MineralIndexMapTest)
{
    // The solver sorts the minerals by name (ANH before CAL), so the deck order CAL, ANH is
    // reversed in the solver. Lower case names are accepted as in the solver.
    const auto interface = make_interface(default_species, {"CAL", "ANH"});

    const std::vector<std::string> deck_minerals = {"CAL", "anh"};
    const auto map = interface->mineral_index_map(deck_minerals);
    BOOST_REQUIRE_EQUAL(map.size(), deck_minerals.size());
    BOOST_CHECK_EQUAL(map[0], 1u);
    BOOST_CHECK_EQUAL(map[1], 0u);

    // Each deck mineral maps to the solver mineral with the same name
    const auto solver_minerals = interface->GetMineralNames();
    BOOST_CHECK_EQUAL(solver_minerals[map[0]], "CAL");
    BOOST_CHECK_EQUAL(solver_minerals[map[1]], "ANH");

    // A mineral not in the solver, and a mineral given twice, are errors
    BOOST_CHECK_EXCEPTION(interface->mineral_index_map({"CAL", "HAL"}),
                          std::runtime_error,
                          message_contains_all({"MINERAL:", "HAL", "not found"}));
    BOOST_CHECK_EXCEPTION(interface->mineral_index_map({"CAL", "cal"}),
                          std::runtime_error,
                          message_contains_all({"MINERAL: duplicate", "cal"}));
}


BOOST_AUTO_TEST_CASE(InitialMineralConcentrationReferenceValuesTest)
{
    const auto interface = make_interface(default_species, {"CAL"});
    const double porosity = 0.25;

    std::vector<double> Cmin_1;
    std::vector<double> Cmin_2;
    interface->calculate_initial_mineral_concentration(Cmin_1, porosity, {{"CAL", 0.001}});
    interface->calculate_initial_mineral_concentration(Cmin_2, porosity, {{"CAL", 0.002}});

    // Not exactly proportional to the weight fraction, since the rock density also depends on it.
    // The reference values depend on the mineral database (molecular weight of CAL, densities).
    BOOST_REQUIRE_EQUAL(Cmin_1.size(), 1);
    BOOST_REQUIRE_EQUAL(Cmin_2.size(), 1);
    BOOST_CHECK_CLOSE(Cmin_1[0], 0.0809281, 1e-3);
    BOOST_CHECK_CLOSE(Cmin_2[0], 0.1618568, 1e-3);
}


BOOST_AUTO_TEST_CASE(InitialMineralConcentrationPorosityDependenceTest)
{
    const auto interface = make_interface(default_species, {"CAL"});

    // The mineral concentration is proportional to (1 - porosity) / porosity, so
    // Cmin * porosity / (1 - porosity) must not depend on the porosity. Comparing against a
    // value computed in the test keeps this independent of the mineral database.
    auto scaled_concentration = [&interface](const double porosity) {
        std::vector<double> Cmin;
        interface->calculate_initial_mineral_concentration(Cmin, porosity, {{"CAL", 0.001}});
        BOOST_REQUIRE_EQUAL(Cmin.size(), 1);
        return Cmin[0] * porosity / (1.0 - porosity);
    };

    const double reference = scaled_concentration(0.25);
    for (const double porosity : {0.1, 0.4}) {
        BOOST_CHECK_CLOSE(scaled_concentration(porosity), reference, 1e-8);
    }

    // Lower porosity gives a higher concentration per pore volume
    std::vector<double> Cmin_low;
    std::vector<double> Cmin_high;
    interface->calculate_initial_mineral_concentration(Cmin_low, 0.1, {{"CAL", 0.001}});
    interface->calculate_initial_mineral_concentration(Cmin_high, 0.4, {{"CAL", 0.001}});
    BOOST_CHECK_GT(Cmin_low[0], Cmin_high[0]);
}


BOOST_AUTO_TEST_CASE(InitialMineralConcentrationResizesOutputTest)
{
    const auto interface = make_interface(default_species, {"CAL"});

    // The output vector is resized to the number of minerals, whatever its previous size
    std::vector<double> Cmin(5, -1.0);
    interface->calculate_initial_mineral_concentration(Cmin, 0.25, {{"CAL", 0.001}});
    BOOST_REQUIRE_EQUAL(Cmin.size(), 1);
    BOOST_CHECK_CLOSE(Cmin[0], 0.0809281, 1e-3);
}


BOOST_AUTO_TEST_CASE(InitialMineralConcentrationUnknownMineralThrowsTest)
{
    const auto interface = make_interface(default_species, {"CAL"});

    std::vector<double> Cmin;
    BOOST_CHECK_EXCEPTION(
        interface->calculate_initial_mineral_concentration(Cmin, 0.25, {{"NOT_A_MINERAL", 0.001}}),
        std::runtime_error,
        message_contains("NOT_A_MINERAL"));

    // Weight fractions for a different number of minerals than the solver has
    BOOST_CHECK_EXCEPTION(interface->calculate_initial_mineral_concentration(Cmin, 0.25, {}),
                          std::runtime_error,
                          message_contains_all({"Weight fractions are given for 0 minerals",
                                                "the geochemistry solver has 1"}));
}


BOOST_AUTO_TEST_CASE(SurfaceConcentrationScalesWithWaterSaturationTest)
{
    const std::vector<std::string> ion_ex = {"X"};
    const auto interface = make_interface(default_species, {"CAL"}, ion_ex);
    const std::size_t nAq = interface->numberOfAqueousBasisSpecies();
    const std::size_t nBasis = interface->numberOfBasisSpecies();
    const double Cion = 1.1e-3;

    // The ion exchanger is inserted in Ctot after the aqueous species, as Cion / swat
    for (const double swat : {1.0, 0.5, 0.25}) {
        std::vector<double> Ctot(nBasis, 0.0);
        double frac_dl = 0.0;
        interface->set_surface_concentrations(swat, Ctot, frac_dl, {{"X", Cion}});
        BOOST_CHECK_CLOSE(Ctot[nAq], Cion / swat, 1e-6);
    }
}


BOOST_AUTO_TEST_CASE(SurfaceConcentrationOnlyChangesIonExchangeEntryTest)
{
    const auto interface = make_interface(default_species, {"CAL"}, {"X"});
    const std::size_t nAq = interface->numberOfAqueousBasisSpecies();

    // Only the ion exchange entry is added to Ctot, the aqueous entries are left as they were
    std::vector<double> Ctot(interface->numberOfBasisSpecies(), 0.0);
    double frac_dl = 0.0;
    interface->set_surface_concentrations(0.5, Ctot, frac_dl, {{"X", 1.1e-3}});
    for (std::size_t i = 0; i < nAq; ++i) {
        BOOST_CHECK_EQUAL(Ctot[i], 0.0);
    }
    BOOST_CHECK_GT(Ctot[nAq], 0.0);
}


BOOST_AUTO_TEST_CASE(SurfaceConcentrationWithoutDiffuseLayerLeavesFractionUnchangedTest)
{
    // The diffusion layer fraction is only assigned when the system has a diffuse layer. With only
    // an ion exchanger it is left as passed in, so a sentinel value must survive the call.
    const auto interface = make_interface(default_species, {"CAL"}, {"X"});

    std::vector<double> Ctot(interface->numberOfBasisSpecies(), 0.0);
    double frac_dl = -1.0;
    interface->set_surface_concentrations(0.5, Ctot, frac_dl, {{"X", 1.1e-3}});
    BOOST_CHECK_EQUAL(frac_dl, -1.0);
}


BOOST_AUTO_TEST_CASE(SurfaceConcentrationUnknownIonExchangeThrowsTest)
{
    const auto interface = make_interface(default_species, {"CAL"}, {"X"});

    std::vector<double> Ctot(interface->numberOfBasisSpecies(), 0.0);
    double frac_dl = 0.0;
    BOOST_CHECK_THROW(
        interface->set_surface_concentrations(0.5, Ctot, frac_dl, {{"NOT_AN_EXCHANGER", 1.0e-3}}),
        std::out_of_range);

    // Concentrations for a different number of ion exchange species than the solver has
    BOOST_CHECK_EXCEPTION(
        interface->set_surface_concentrations(0.5, Ctot, frac_dl, {}),
        std::runtime_error,
        message_contains_all({"Concentrations are given for 0 ion exchange species",
                              "the geochemistry solver has 1"}));
}

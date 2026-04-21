#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <fstream>
#include <sstream>

#include <opm/simulators/geochemistry/StandaloneSolvers.hpp>
#include <opm/simulators/geochemistry/Thermo/hkf.h>
#include <opm/simulators/geochemistry/Thermo/eps_JN.h>
#include <opm/simulators/geochemistry/Thermo/ions.h>
#include <opm/simulators/geochemistry/Thermo/ThermoTable.h>
#include <opm/simulators/geochemistry/Thermo/thermodata.h>
#include <opm/simulators/geochemistry/Thermo/water.h>

static constexpr double abs_tolerance_ = 1.0e-7;
static constexpr double abs_tolerance2_ = 1.0e-4;  // for most water props.

/*
* Stores reference cases 1, 2, and 3 from Table 5 of the IAPWS-97
* industrial thermodynamic formulation.
*
* See also:
*
*   - https://github.com/jjgomera/iapws.git
*   - http://twt.mpei.ac.ru/mcs/worksheets/iapws/IAPWS-IF97-Region1.xmcd
*   - http://www.casprod.eu/~leon/research/iapws-if97/listing.html
*/
struct WaterPropertiesIAPWS_97
{
    WaterPropertiesIAPWS_97(int reference_case)
    {
        switch (reference_case) {
            case 1:
            {
                pressure_ = 3.0e6;
                temperature_ = 300.0;
                v_ = 0.100215168e-2;
                h_ = 0.115331273e3;
                u_ = 0.112324818e3;
                s_ = 0.392294792;
                cp_ = 0.417301218e1;
                w_ = 0.150773921e4;
                break;
            }
            case 2:
            {
                pressure_ = 80.0e6;
                temperature_ = 300.0;
                v_ = 0.971180894e-3;
                h_ = 0.184142828e3;
                u_ = 0.106448356e3;
                s_ = 0.368563852;
                cp_ = 0.401008987e1;
                w_ = 0.163469054e4;
                break;
            }
            case 3:
            {
                pressure_ = 3.0e6;
                temperature_ = 500.0;
                v_ = 0.120241800e-2;
                h_ = 0.975542239e3;
                u_ = 0.971934985e3;
                s_ = 0.258041912e1;
                cp_ = 0.465580682e1;
                w_ = 0.124071337e4;
                break;
            }
        }
    }

    double pressure_;
    double temperature_;

    double v_;  // specific volume (m3/kg)
    double h_;  // specific enthalpy (kJ/kg)
    double u_;  // specific internal energy (kJ/kg)
    double s_;  // specific entropy (kJ/kg/K)
    double cp_;  // specific isobaric heat capacity (kJ/kg/K)
    double w_;  // speed of sound (m/s)

    //double cv_;  // specific isochoric heat capacity (kJ/kg/K)
};


TEST_CASE("Test calculation of thermodynamic properties for water")
{

    water water_props;

    static const std::array<WaterPropertiesIAPWS_97, 3> reference_cases =
    {
        WaterPropertiesIAPWS_97(1),
        WaterPropertiesIAPWS_97(2),
        WaterPropertiesIAPWS_97(3),
    };

    for(const auto& ref_props: reference_cases)
    {
        const double T = ref_props.temperature_;
        const double P = ref_props.pressure_;

        water_props.gibbsIAPWS(T, P);

        CHECK_THAT(water_props.v_, Catch::Matchers::WithinAbs(ref_props.v_, abs_tolerance_));

        // Note: Reference results are in kJ/kg or kJ/kg/K, we use Joules...
        CHECK_THAT(1.0e-3*water_props.h_, Catch::Matchers::WithinAbs(ref_props.h_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.u_, Catch::Matchers::WithinAbs(ref_props.u_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.s_, Catch::Matchers::WithinAbs(ref_props.s_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.cp_, Catch::Matchers::WithinAbs(ref_props.cp_, abs_tolerance2_));
        CHECK_THAT(water_props.w_, Catch::Matchers::WithinAbs(ref_props.w_, abs_tolerance2_));

        water_props.printProperties();
    }

}

TEST_CASE("Test calculation of saturation pressure for water")
{
    water water_props;

    // Table 35 of IAPWS97 paper
    static constexpr std::array<std::pair<double, double>, 3> table35_tests =
    {
        std::make_pair<double, double>(300.0, 0.353658941e-2),
        std::make_pair<double, double>(500.0, 0.263889776e1),
        std::make_pair<double, double>(600.0, 0.123443146e2)
    };

    for(const auto& [T, Psat] : table35_tests)
    {
        const double satPress = 1.0e-6*water_props.PsatIAPWS(T);
        CHECK_THAT(satPress, Catch::Matchers::WithinAbs(Psat, abs_tolerance_));
    }
}

TEST_CASE("Test HKF standard-state water properties")
{
    hkf hkf_props;
    hkf_props.WaterProp(298.15, 1.0e5);

    CHECK_THAT(hkf_props.rhow_, Catch::Matchers::WithinAbs(997.047435408, 1.0e-6));
    CHECK_THAT(hkf_props.epsw_, Catch::Matchers::WithinAbs(78.2438908377, 1.0e-6));
    CHECK_THAT(hkf_props.G_, Catch::Matchers::WithinAbs(-237181.719094, 1.0e-6));
    CHECK_THAT(hkf_props.H_, Catch::Matchers::WithinAbs(-285830.819484, 1.0e-6));
    CHECK_THAT(hkf_props.S_, Catch::Matchers::WithinAbs(69.9280637281, 1.0e-6));
    CHECK_THAT(hkf_props.V_, Catch::Matchers::WithinAbs(1.80686287936e-05, 1.0e-12));
    CHECK_THAT(hkf_props.Cp_, Catch::Matchers::WithinAbs(75.3381005398, 1.0e-6));
}

TEST_CASE("Thermo table calculator matches water table-style output")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("H2O", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.T, Catch::Matchers::WithinAbs(25.0, abs_tolerance_));
    CHECK_THAT(row.P, Catch::Matchers::WithinAbs(1.0, abs_tolerance_));
    CHECK_THAT(row.rho, Catch::Matchers::WithinAbs(0.9970614, 2.0e-5));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(41.55238, 1.0e-3));
    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-237181.4, 1.0e1));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-285837.3, 1.0e1));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(69.92418, 1.0e-2));
    CHECK_THAT(row.V, Catch::Matchers::WithinAbs(18.06830, 1.0e-2));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(75.36053, 5.0e-2));
}

TEST_CASE("Thermo table calculator broadcasts vector inputs")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("H2O", {25.0, 50.0}, {1.0});

    REQUIRE(rows.size() == 2);
    CHECK_THAT(rows[0].T, Catch::Matchers::WithinAbs(25.0, abs_tolerance_));
    CHECK_THAT(rows[1].T, Catch::Matchers::WithinAbs(50.0, abs_tolerance_));
    CHECK_THAT(rows[0].P, Catch::Matchers::WithinAbs(1.0, abs_tolerance_));
    CHECK_THAT(rows[1].P, Catch::Matchers::WithinAbs(1.0, abs_tolerance_));
    CHECK(rows[1].rho < rows[0].rho);
}

TEST_CASE("Thermo table calculator formats table output")
{
    ThermoTableCalculator calculator;
    const auto table = calculator.evaluateFormatted("H2O", {25.0, 50.0}, {1.0});

    CHECK(table.find("T") != std::string::npos);
    CHECK(table.find("P") != std::string::npos);
    CHECK(table.find("rho") != std::string::npos);
    CHECK(table.find("Cp") != std::string::npos);
    CHECK(table.find("   1   25.00") != std::string::npos);
    CHECK(table.find("   2   50.00") != std::string::npos);
    CHECK(table.find("75.33810") != std::string::npos);
    CHECK(table.find("75.29589") != std::string::npos);
}

TEST_CASE("Thermo table calculator exposes direct IAPWS water saturation pressure")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluateWaterSaturationPressure({26.85, 226.85});

    REQUIRE(rows.size() == 2);
    CHECK_THAT(rows[0].T, Catch::Matchers::WithinAbs(26.85, abs_tolerance_));
    CHECK_THAT(rows[0].Psat, Catch::Matchers::WithinAbs(0.0353658941, 1.0e-8));
    CHECK_THAT(rows[1].Psat, Catch::Matchers::WithinAbs(26.3889776, 1.0e-6));

    const auto formatted = calculator.evaluateWaterSaturationPressureFormatted({25.0});
    CHECK(formatted.find("Psat") != std::string::npos);
}

TEST_CASE("Thermo table calculator reproduces HKF basis species reference data")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("Na+", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-62591.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-57433.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(13.96*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.V, Catch::Matchers::WithinAbs(-1206.704398355, 1.0e-6));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(38.1192278706, 1.0e-4));
}

TEST_CASE("Thermo table calculator reproduces HKF mineral reference data")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("CALCITE", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-269880.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-288552.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(22.15*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.V, Catch::Matchers::WithinAbs(36.934, abs_tolerance2_));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(81.876, 2.0e-1));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(1.84865, 1.0e-4));
}

TEST_CASE("Thermo table calculator reproduces HKF gas reference data")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("CO2,g", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    CHECK_THAT(row.G, Catch::Matchers::WithinAbs(-94254.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.H, Catch::Matchers::WithinAbs(-94051.0*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.S, Catch::Matchers::WithinAbs(51.085*UnitConversionFactors::cal2J_, abs_tolerance2_));
    CHECK_THAT(row.Cp, Catch::Matchers::WithinAbs(37.155, 2.0e-1));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(-7.81373, 1.0e-4));
}

TEST_CASE("Thermo table calculator derives H2O,g logK from IAPWS saturation pressure")
{
    ThermoTableCalculator calculator;
    const auto rows = calculator.evaluate("H2O,g", {25.0}, {1.0});

    REQUIRE(rows.size() == 1);
    const auto& row = rows.front();

    const double expected_logK = std::log10(PhysicalConstants::standard_gas_pressure / water::PsatIAPWS(298.15));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(expected_logK, 1.0e-10));
    CHECK_THAT(row.logK, Catch::Matchers::WithinAbs(1.49897, 1.0e-4));
}

TEST_CASE("Thermo table standalone solver reads debug input")
{
    const std::string input = R"(SPECIESLIST
H2O
CALCITE
/ end
TEMPS
25.0 50.0
/ end
PRESSURES
1.0e5
/ end
TempUnit C
PresUnit Pa
IncludeIndex 1
)";

    std::istringstream input_stream(input);
    ThermoTableSolver solver;
    const auto output = solver.solve("thermotable_debug", input_stream);

    CHECK(output.find("SPECIES H2O") != std::string::npos);
    CHECK(output.find("SPECIES CALCITE") != std::string::npos);
    CHECK(output.find("IAPWS97_PSAT_H2O") != std::string::npos);
    CHECK(output.find("   1   25.00") != std::string::npos);
    CHECK(output.find("   2   50.00") != std::string::npos);
    CHECK(output.find("18.06863") != std::string::npos);
    CHECK(output.find("36.93400") != std::string::npos);
}

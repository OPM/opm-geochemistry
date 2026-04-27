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
#ifndef THERMO_TABLE_H
#define THERMO_TABLE_H

#include <string>
#include <vector>

#include <opm/simulators/geochemistry/Core/ChemTable.h>
#include <opm/simulators/geochemistry/Thermo/hkf.h>

struct ThermoTableRow
{
    double T{0.0};     // Celsius
    double P{0.0};     // bar
    double rho{0.0};   // g / cm^3
    double logK{0.0};  // dimensionless
    double G{0.0};     // J / mol
    double H{0.0};     // J / mol
    double S{0.0};     // J / mol / K
    double V{0.0};     // cm^3 / mol
    double Cp{0.0};    // J / mol / K
};

struct WaterSaturationPressureRow
{
    double T{0.0};     // Celsius
    double Psat{0.0};  // bar
};

class ThermoTableCalculator
{
public:
    ThermoTableCalculator();

    [[nodiscard]] std::vector<ThermoTableRow> evaluate(const std::string& species_name,
                                                       const std::vector<double>& temperatures_celsius,
                                                       const std::vector<double>& pressures_bar = {});
    [[nodiscard]] std::string evaluateFormatted(const std::string& species_name,
                                                const std::vector<double>& temperatures_celsius,
                                                const std::vector<double>& pressures_bar = {},
                                                bool include_index = true);
    [[nodiscard]] static std::string formatRows(const std::vector<ThermoTableRow>& rows,
                                                bool include_index = true);
    [[nodiscard]] std::vector<WaterSaturationPressureRow> evaluateWaterSaturationPressure(const std::vector<double>& temperatures_celsius);
    [[nodiscard]] std::string evaluateWaterSaturationPressureFormatted(const std::vector<double>& temperatures_celsius,
                                                                       bool include_index = true);
    [[nodiscard]] static std::string formatWaterSaturationPressureRows(const std::vector<WaterSaturationPressureRow>& rows,
                                                                       bool include_index = true);
    [[nodiscard]] static bool isWaterSpeciesName(const std::string& species_name);

private:
    enum class SpeciesKind
    {
        Water,
        Basis,
        Complex,
        Mineral
    };

    struct SpeciesRef
    {
        SpeciesKind kind{SpeciesKind::Basis};
        int index{-1};
    };

    ChemTable basis_;
    ChemTable aqueous_;
    ChemTable minerals_;
    hkf eos_;

    [[nodiscard]] SpeciesRef findSpecies(const std::string& species_name) const;
    [[nodiscard]] static int findRowIndex(const ChemTable& table, const std::string& species_name);
    [[nodiscard]] static std::vector<double> broadcastInput(const std::vector<double>& values,
                                                            std::size_t size,
                                                            double default_value);
    [[nodiscard]] ThermoTableRow evaluateOne(const SpeciesRef& species, double temperature_celsius, double pressure_bar);
    [[nodiscard]] double reactionLogK(const ChemTable& table,
                                      int row_index,
                                      const StandardStateProperties& species_props,
                                      double temperature_kelvin,
                                      double pressure_pa);
    [[nodiscard]] StandardStateProperties basisProperties(int row_index, double temperature_kelvin, double pressure_pa);
};

#endif

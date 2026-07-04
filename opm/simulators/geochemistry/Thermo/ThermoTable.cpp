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
#include <opm/simulators/geochemistry/Thermo/ThermoTable.h>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

constexpr double pa_per_bar = 1.0e5;
constexpr double bar_per_pa = 1.0e-5;
constexpr double cm3_per_m3 = 1.0e6;
constexpr double g_cm3_per_kg_m3 = 1.0e-3;

void appendField(std::ostringstream& stream, double value, int width, int precision)
{
    stream << std::setw(width) << std::fixed << std::setprecision(precision) << value;
}

} // namespace

ThermoTableCalculator::ThermoTableCalculator()
{
    basis_.initialize_full_database_hkf(ChemTable::Type::Basis);
    basis_.init_tables();
    basis_.update_number_of_analytical_species();
    basis_.fix_hkf_units(ChemTable::Type::Basis);

    aqueous_.initialize_full_database_hkf(ChemTable::Type::Complex);
    aqueous_.init_tables();
    aqueous_.update_number_of_analytical_species();
    aqueous_.fix_hkf_units(ChemTable::Type::Complex);

    minerals_.initialize_full_database_hkf(ChemTable::Type::Mineral);
    minerals_.init_tables();
    minerals_.update_number_of_analytical_species();
    minerals_.fix_hkf_units(ChemTable::Type::Mineral);
}

std::vector<ThermoTableRow> ThermoTableCalculator::evaluate(const std::string& species_name,
                                                            const std::vector<double>& temperatures_celsius,
                                                            const std::vector<double>& pressures_bar)
{
    if (temperatures_celsius.empty() && pressures_bar.empty())
    {
        return {};
    }

    const std::size_t size = std::max(temperatures_celsius.size(), pressures_bar.size());
    if (size == 0)
    {
        return {};
    }

    const auto temperatures = broadcastInput(temperatures_celsius, size, 25.0);
    const auto pressures = broadcastInput(pressures_bar, size, 1.0);
    const auto species = findSpecies(species_name);

    std::vector<ThermoTableRow> rows;
    rows.reserve(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        rows.push_back(evaluateOne(species, temperatures[i], pressures[i]));
    }
    return rows;
}

std::string ThermoTableCalculator::evaluateFormatted(const std::string& species_name,
                                                     const std::vector<double>& temperatures_celsius,
                                                     const std::vector<double>& pressures_bar,
                                                     bool include_index)
{
    return formatRows(evaluate(species_name, temperatures_celsius, pressures_bar), include_index);
}

std::string ThermoTableCalculator::formatRows(const std::vector<ThermoTableRow>& rows, bool include_index)
{
    if (rows.empty())
    {
        return {};
    }

    std::ostringstream stream;
    if (include_index)
    {
        stream << std::setw(4) << "";
    }

    stream << std::setw(8) << "T"
           << std::setw(11) << "P"
           << std::setw(10) << "rho"
           << std::setw(10) << "logK"
           << std::setw(11) << "G"
           << std::setw(11) << "H"
           << std::setw(10) << "S"
           << std::setw(10) << "V"
           << std::setw(10) << "Cp"
           << '\n';

    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        const auto& row = rows[i];
        if (include_index)
        {
            stream << std::setw(4) << (i + 1);
        }

        appendField(stream, row.T, 8, 2);
        appendField(stream, row.P, 11, 6);
        appendField(stream, row.rho, 10, 7);
        appendField(stream, row.logK, 10, 5);
        appendField(stream, row.G, 11, 1);
        appendField(stream, row.H, 11, 1);
        appendField(stream, row.S, 10, 5);
        appendField(stream, row.V, 10, 5);
        appendField(stream, row.Cp, 10, 5);
        stream << '\n';
    }

    return stream.str();
}

std::vector<WaterSaturationPressureRow> ThermoTableCalculator::evaluateWaterSaturationPressure(const std::vector<double>& temperatures_celsius)
{
    if (temperatures_celsius.empty())
    {
        return {};
    }

    std::vector<WaterSaturationPressureRow> rows;
    rows.reserve(temperatures_celsius.size());
    for (const double temperature_celsius : temperatures_celsius)
    {
        WaterSaturationPressureRow row;
        row.T = temperature_celsius;
        row.Psat = eos_.waterSaturationPressure(temperature_celsius + 273.15)*bar_per_pa;
        rows.push_back(row);
    }
    return rows;
}

std::string ThermoTableCalculator::evaluateWaterSaturationPressureFormatted(const std::vector<double>& temperatures_celsius,
                                                                            bool include_index)
{
    return formatWaterSaturationPressureRows(evaluateWaterSaturationPressure(temperatures_celsius), include_index);
}

std::string ThermoTableCalculator::formatWaterSaturationPressureRows(const std::vector<WaterSaturationPressureRow>& rows,
                                                                     bool include_index)
{
    if (rows.empty())
    {
        return {};
    }

    std::ostringstream stream;
    if (include_index)
    {
        stream << std::setw(4) << "";
    }

    stream << std::setw(8) << "T"
           << std::setw(12) << "Psat"
           << '\n';

    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        const auto& row = rows[i];
        if (include_index)
        {
            stream << std::setw(4) << (i + 1);
        }

        appendField(stream, row.T, 8, 2);
        appendField(stream, row.Psat, 12, 6);
        stream << '\n';
    }

    return stream.str();
}

bool ThermoTableCalculator::isWaterSpeciesName(const std::string& species_name)
{
    return species_name == "H2O"
        || species_name == "water"
        || species_name == "WATER"
        || species_name == "H2O,g"
        || species_name == "H2O(G)";
}

ThermoTableCalculator::SpeciesRef ThermoTableCalculator::findSpecies(const std::string& species_name) const
{
    if (species_name == "H2O" || species_name == "water" || species_name == "WATER")
    {
        return {SpeciesKind::Water, basis_.get_row_index("H2O")};
    }

    if (const int index = findRowIndex(basis_, species_name); index >= 0)
    {
        return {SpeciesKind::Basis, index};
    }
    if (const int index = findRowIndex(aqueous_, species_name); index >= 0)
    {
        return {SpeciesKind::Complex, index};
    }
    if (const int index = findRowIndex(minerals_, species_name); index >= 0)
    {
        return {SpeciesKind::Mineral, index};
    }

    throw std::invalid_argument("Unknown species name: " + species_name);
}

int ThermoTableCalculator::findRowIndex(const ChemTable& table, const std::string& species_name)
{
    int index = table.get_row_index(species_name);
    if (index >= 0)
    {
        return index;
    }

    const auto it = std::find(table.nick_name_.begin(), table.nick_name_.end(), species_name);
    if (it != table.nick_name_.end())
    {
        return static_cast<int>(std::distance(table.nick_name_.begin(), it));
    }
    return -1;
}

std::vector<double> ThermoTableCalculator::broadcastInput(const std::vector<double>& values,
                                                          std::size_t size,
                                                          double default_value)
{
    if (values.empty())
    {
        return std::vector<double>(size, default_value);
    }
    if (values.size() == size)
    {
        return values;
    }
    if (values.size() == 1)
    {
        return std::vector<double>(size, values.front());
    }

    throw std::invalid_argument("Temperature and pressure vectors must either have matching sizes or size 1.");
}

ThermoTableRow ThermoTableCalculator::evaluateOne(const SpeciesRef& species, double temperature_celsius, double pressure_bar)
{
    const double temperature_kelvin = temperature_celsius + 273.15;
    const double pressure_pa = pressure_bar*pa_per_bar;

    StandardStateProperties props;
    double logK = 0.0;

    if (species.kind == SpeciesKind::Water)
    {
        props = eos_.waterProperties(temperature_kelvin, pressure_pa);
        logK = -props.G / (NumericalConstants::LNTEN*PhysicalConstants::IdealGasConstant*temperature_kelvin);
    }
    else if (species.kind == SpeciesKind::Basis)
    {
        props = basisProperties(species.index, temperature_kelvin, pressure_pa);
        logK = -props.G / (NumericalConstants::LNTEN*PhysicalConstants::IdealGasConstant*temperature_kelvin);
    }
    else if (species.kind == SpeciesKind::Complex)
    {
        if (aqueous_.model_[species.index] != LogKModel::HKF)
        {
            throw std::invalid_argument("Only HKF aqueous species are currently supported in ThermoTableCalculator.");
        }
        props = eos_.ionProperties(temperature_kelvin,
                                   pressure_pa,
                                   aqueous_.deltaG_[species.index],
                                   aqueous_.deltaH_[species.index],
                                   aqueous_.S_[species.index],
                                   aqueous_.a1_[species.index],
                                   aqueous_.a2_[species.index],
                                   aqueous_.a3_[species.index],
                                   aqueous_.a4_[species.index],
                                   aqueous_.c1_[species.index],
                                   aqueous_.c2_[species.index],
                                   aqueous_.omega_[species.index],
                                   aqueous_.charge_[species.index],
                                   aqueous_.re_[species.index]);
        logK = reactionLogK(aqueous_, species.index, props, temperature_kelvin, pressure_pa);
    }
    else
    {
        if (minerals_.model_[species.index] != LogKModel::HKF)
        {
            throw std::invalid_argument("Only HKF minerals/gases are currently supported in ThermoTableCalculator.");
        }
        props = eos_.mineralProperties(temperature_kelvin,
                                       pressure_pa,
                                       minerals_.deltaG_[species.index],
                                       minerals_.deltaH_[species.index],
                                       minerals_.S_[species.index],
                                       minerals_.a1_[species.index],
                                       minerals_.a2_[species.index],
                                       minerals_.a3_[species.index],
                                       minerals_.mol_volume_[species.index]);
        logK = reactionLogK(minerals_, species.index, props, temperature_kelvin, pressure_pa);
        if (minerals_.row_name_[species.index] == "H2O,g")
        {
            logK = eos_.waterVaporLogK(temperature_kelvin, PhysicalConstants::standard_gas_pressure);
        }
    }

    ThermoTableRow row;
    row.T = temperature_celsius;
    row.P = pressure_bar;
    row.rho = eos_.rhow_*g_cm3_per_kg_m3;
    row.logK = logK;
    row.G = props.G;
    row.H = props.H;
    row.S = props.S;
    row.V = props.V*cm3_per_m3;
    row.Cp = props.Cp;
    return row;
}

double ThermoTableCalculator::reactionLogK(const ChemTable& table,
                                           int row_index,
                                           const StandardStateProperties& species_props,
                                           double temperature_kelvin,
                                           double pressure_pa)
{
    double logK = species_props.G;
    for (int column_index = 0; column_index < table.noColumns_; ++column_index)
    {
        const double stoichiometric_coefficient = table.M_[row_index][column_index];
        if (stoichiometric_coefficient == 0.0)
        {
            continue;
        }

        const int basis_index = basis_.get_row_index(table.col_name_[column_index]);
        if (basis_index < 0)
        {
            continue;
        }

        const auto basis_props = basisProperties(basis_index, temperature_kelvin, pressure_pa);
        logK -= stoichiometric_coefficient*basis_props.G;
    }

    logK /= NumericalConstants::LNTEN*PhysicalConstants::IdealGasConstant*temperature_kelvin;
    return logK;
}

StandardStateProperties ThermoTableCalculator::basisProperties(int row_index, double temperature_kelvin, double pressure_pa)
{
    if (basis_.row_name_[row_index] == "H2O")
    {
        return eos_.waterProperties(temperature_kelvin, pressure_pa);
    }
    if (basis_.model_[row_index] != LogKModel::HKF)
    {
        throw std::invalid_argument("Only HKF basis species are currently supported in ThermoTableCalculator.");
    }

    return eos_.ionProperties(temperature_kelvin,
                              pressure_pa,
                              basis_.deltaG_[row_index],
                              basis_.deltaH_[row_index],
                              basis_.S_[row_index],
                              basis_.a1_[row_index],
                              basis_.a2_[row_index],
                              basis_.a3_[row_index],
                              basis_.a4_[row_index],
                              basis_.c1_[row_index],
                              basis_.c2_[row_index],
                              basis_.omega_[row_index],
                              basis_.charge_[row_index],
                              basis_.re_[row_index]);
}

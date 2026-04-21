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
/**
 * @file
 */
#ifndef CONSTANTS_IS_DEF_HPP
#define CONSTANTS_IS_DEF_HPP

struct PhysicalConstants
{
    static constexpr double Faraday{9.648456e+4};  // C / mol
    static constexpr double VacuumPermittivity{8.854187817e-12};  // A s V^-1 m^-1
    static constexpr double ElectronCharge{1.6021766208e-19};  // C
    static constexpr double Avogadro{6.0221415e+23};  // 1/mol
    static constexpr double Boltzmann{1.3806503e-23};  // m^2 kg s^-2 K^-1
    static constexpr double IdealGasConstant{8.31441}; //  J / K / mol (R_g = N_avo*k_Boltzmann)

    static constexpr double atmospheric_pressure{1.01325e5}; // Pa
    static constexpr double standard_gas_pressure{1.0e5}; // Pa, gas activity reference pressure
    static constexpr double ambient_temperature{ 298.15 }; // K

};

struct UnitConversionFactors
{
    static constexpr double cal2J_{4.184};
    static constexpr double bar2Pa_{1e5};
    static constexpr double eta_{1.66027}; //  ?cal/mol*10^5
};

struct NumericalConstants
{
    static constexpr double LNTEN{2.302585092994};
    static constexpr double ZERO_THRESHOLD{1.0e-16};
    static constexpr double ZERO_CONCENTRATION_THRESHOLD{1e-9}; // Used by 1D transport solver.
    static constexpr double LOG10_ALMOST_ZERO{-20.0};

    static constexpr double CONCENTRATION_THRESHOLD = 10.0; // 10 mol/L
    static constexpr double POTENTIAL_THRESHOLD = 1.0;  // potential of 1 V
};

struct DefaultValues
{
    static constexpr double MolarConcentration{1.0e-8};
    static constexpr double pH{7.0};
    static constexpr double pe{ 4.0 };
    static constexpr double Temperature{298.15};  // 25 degrees Celsius
    static constexpr double Pressure{1.0e6};  // 10 bar
    static constexpr double WaterDielectricConstant{80.0};
    static constexpr double WaterDensity{1000.0};  // Unit: [kg/m3]
};

#endif  // CONSTANTS_IS_DEF_HPP

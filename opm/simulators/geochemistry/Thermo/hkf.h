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
#ifndef HKF_H
#define HKF_H

#include <memory>

#include <opm/simulators/geochemistry/Thermo/eps_JN.h>
#include <opm/simulators/geochemistry/Thermo/ions.h>
#include <opm/simulators/geochemistry/Thermo/thermodata.h>
#include <opm/simulators/geochemistry/Thermo/water.h>

struct StandardStateProperties
{
    double G{0.0};   // Gibbs free energy [J/mol]
    double H{0.0};   // Enthalpy [J/mol]
    double S{0.0};   // Entropy [J/mol/K]
    double V{0.0};   // Molar volume [m^3/mol]
    double Cp{0.0};  // Isobaric heat capacity [J/mol/K]
};

class hkf
{
public:
    
    hkf();
    
    void WaterProp(double T, double P);
    [[nodiscard]] double waterSaturationPressure(double T) const;
    [[nodiscard]] double waterVaporLogK(double T, double gas_activity_reference_pressure) const;

    [[nodiscard]] StandardStateProperties waterProperties(double T, double P);
    [[nodiscard]] StandardStateProperties mineralProperties(double T,
                                                            double P,
                                                            double G,
                                                            double H,
                                                            double S,
                                                            double a,
                                                            double b,
                                                            double c,
                                                            double V);
    [[nodiscard]] StandardStateProperties ionProperties(double T,
                                                        double P,
                                                        double G,
                                                        double H,
                                                        double S,
                                                        double a1,
                                                        double a2,
                                                        double a3,
                                                        double a4,
                                                        double c1,
                                                        double c2,
                                                        double omega,
                                                        double Z,
                                                        double re_ref);
    
    void dGMineral(double T, double P, double* G, double* S, double* a, double* b, double* c, double* V, int size, double* logK);
    
    void dGIons(double T, double P, double* G, double* S,
                double* a1, double* a2, double* a3, double* a4,
                double* c1, double* c2,
                double* omega, double* Z, double* re_ref, int size,
                int skip, double* dG, double* MV);
    
    double rhow_;
    double epsw_;
    double G_;   // Standard molal Gibbs free energy of formation for H2O [J/mol]
    double H_;   // Standard molal enthalpy of formation for H2O [J/mol]
    double S_;   // Standard molal entropy for H2O [J/mol/K]
    double V_;   // Standard molal volume for H2O [m^3/mol]
    double Cp_;  // Standard molal isobaric heat capacity for H2O [J/mol/K]
    
private:
    
    std::unique_ptr<water> W_;
    std::unique_ptr<eps_JN> EPS_;
    std::unique_ptr<ions> BORN_;
    
    static constexpr double Rg_inv_ = 0.12027239580856473682150703284415;
    static constexpr double loge_ = 0.43429448190325182765112891891661;
    static constexpr double ln10_ = 2.3025850929940456840179914546844;
    static constexpr double water_molar_mass_ = 18.01528e-3;  // kg/mol
    
    double Tref_;
    double Tref2_inv_;
    double Pref_;
    double Pref_inv_;
    double Theta_;
    double Psi_;
    
    double Y_ref_;
    double Z_ref_;
    
    double Gtr_;
    double Htr_;
    double Str_;
    double Ttr_;

    void updateWaterState(double T, double P);
    [[nodiscard]] static double propertyShift(double G_ref, double H_ref, double S_ref, double T_ref);
    
};

#endif

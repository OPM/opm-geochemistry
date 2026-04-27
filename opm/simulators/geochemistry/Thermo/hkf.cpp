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
#include <opm/simulators/geochemistry/Thermo/hkf.h>

#include <cmath>

hkf::hkf()
: rhow_(0.0)
, epsw_(0.0)
, G_(0.0)
, H_(0.0)
, S_(0.0)
, V_(0.0)
, Cp_(0.0)
, Tref_(298.15)  // K
, Tref2_inv_(1.0/88893.4225) // K^-2
, Pref_(1.0e5)  // Pa
, Pref_inv_(1.0e-5) // Pa^-1
, Theta_(228.0)  // K
, Psi_(2600.0e5)  // Pa
, Y_ref_(-5.7950892275586268e-005)
, Z_ref_(-0.012780550523412457)
, Gtr_(-235517.360)  // J / mol
, Htr_(-287721.128) // J / mol
, Str_(63.312288)  // J / mol / K
, Ttr_(273.16)
{
    W_ = std::make_unique<water>();
    EPS_ = std::make_unique<eps_JN>(W_.get());
    BORN_ = std::make_unique<ions>(W_.get());
    
    /*
     // Calculated from:
     W_->gibbsIAPWS(Tref_, Pref_);
     EPS_->permittivity_TP(Tref_, Pref_);
     Y_ref_= EPS_->bornY_;
     BORN_->born_df(Tref_, Pref_);
     */
}

double hkf::propertyShift(double G_ref, double H_ref, double S_ref, double T_ref)
{
    return H_ref - G_ref - T_ref*S_ref;
}

void hkf::updateWaterState(double T, double P)
{
    EPS_->permittivity_TP(T, P);
    rhow_ = W_->denst_;
    epsw_ = EPS_->permittivity_;
}

void hkf::WaterProp(double T, double P)
{
    const auto props = waterProperties(T, P);
    G_ = props.G;
    H_ = props.H;
    S_ = props.S;
    V_ = props.V;
    Cp_ = props.Cp;
}

double hkf::waterSaturationPressure(double T) const
{
    return water::PsatIAPWS(T);
}

double hkf::waterVaporLogK(double T, double gas_activity_reference_pressure) const
{
    return std::log10(gas_activity_reference_pressure / waterSaturationPressure(T));
}

StandardStateProperties hkf::waterProperties(double T, double P)
{
    updateWaterState(T, P);

    StandardStateProperties props;
    props.G = W_->G_ + Gtr_ + Ttr_*Str_ - T*Str_;
    props.H = W_->H_ + Htr_;
    props.S = W_->s_ * water_molar_mass_ + Str_;
    props.V = W_->v_ * water_molar_mass_;
    props.Cp = W_->cp_ * water_molar_mass_;
    return props;
}

StandardStateProperties hkf::mineralProperties(double T,
                                               double P,
                                               double G_ref,
                                               double H_ref,
                                               double S_ref,
                                               double a,
                                               double b,
                                               double c,
                                               double V_ref)
{
    updateWaterState(T, P);

    const double dT = T - Tref_;
    const double dP = P - Pref_;
    const double T_inv = 1.0 / T;
    const double T_inv2 = T_inv*T_inv;
    const double logTTref = std::log(T / Tref_);

    StandardStateProperties props;
    props.G = G_ref - S_ref*dT
              + a*(dT - T*logTTref)
              - 0.5*b*dT*dT
              - 0.5*c*Tref2_inv_*T_inv*dT*dT
              + V_ref*dP;
    props.S = S_ref + a*logTTref + b*dT + 0.5*c*(Tref2_inv_ - T_inv2);
    props.Cp = a + b*T + c*T_inv2;
    props.V = V_ref;
    props.H = props.G + T*props.S + propertyShift(G_ref, H_ref, S_ref, Tref_);
    return props;
}

StandardStateProperties hkf::ionProperties(double T,
                                           double P,
                                           double G_ref,
                                           double H_ref,
                                           double S_ref,
                                           double a1,
                                           double a2,
                                           double a3,
                                           double a4,
                                           double c1,
                                           double c2,
                                           double omega,
                                           double Z,
                                           double re_ref)
{
    updateWaterState(T, P);
    BORN_->born_df(T, P);

    const double dT = T - Tref_;
    const double Tp = T / Tref_;
    const double dP = P*Pref_inv_ - 1.0;

    const double theta_diff = T - Theta_;
    const double theta_diff_inv = 1.0 / theta_diff;
    const double theta_diff_inv2 = theta_diff_inv*theta_diff_inv;
    const double theta_diff_inv3 = theta_diff_inv2*theta_diff_inv;

    const double C1 = dT - T*std::log(T / Tref_);
    const double log_kernel = std::log(Tp*(Tref_ - Theta_) / theta_diff);
    const double C2 = -T / Theta_ / Theta_*log_kernel - dT / Theta_ / (Tref_ - Theta_);
    const double A1 = dP;
    const double A2 = std::log((Psi_ + P) / (Psi_ + Pref_));
    const double A3 = dP*theta_diff_inv;
    const double A4 = A2*theta_diff_inv;

    const double Wref = Z_ref_ - EPS_->bornZ_ + Y_ref_*dT;
    const double W = -EPS_->bornZ_ - 1.0;

    double Wi = 0.0;
    double w_T = 0.0;
    double w_TT = 0.0;
    double w_P = 0.0;
    if (Z != 0.0)
    {
        BORN_->born(Z, re_ref, Wi, w_T, w_TT, w_P);
    }

    StandardStateProperties props;
    props.G = G_ref - S_ref*dT + c1*C1 + c2*C2 + a1*A1 + a2*A2 + a3*A3 + a4*A4 + (Wi - omega)*W + omega*Wref;

    const double C2_T = -log_kernel / (Theta_*Theta_) + 1.0 / (Theta_*theta_diff) - 1.0 / (Theta_*(Tref_ - Theta_));
    props.S = S_ref
              + c1*std::log(T / Tref_)
              - c2*C2_T
              + a3*dP*theta_diff_inv2
              + a4*A2*theta_diff_inv2
              - w_T*W
              + Wi*EPS_->bornY_
              - omega*Y_ref_;

    const double born_TT = w_TT*W - 2.0*w_T*EPS_->bornY_ - Wi*EPS_->bornX_;
    props.Cp = c1
               + c2*theta_diff_inv2
               - 2.0*T*(a3*dP + a4*A2)*theta_diff_inv3
               - T*born_TT;

    const double MV1 = 1.0 / (Psi_*Pref_inv_ + P*Pref_inv_);
    const double Chat = 41.84e-3 / UnitConversionFactors::cal2J_;
    props.V = a1 + a2*MV1 + (a3 + a4*MV1)*theta_diff_inv - 1e5*omega*EPS_->bornQ_ - (EPS_->bornZ_ + 1.0)*w_P*1e5;
    props.V *= Chat;

    props.H = props.G + T*props.S + propertyShift(G_ref, H_ref, S_ref, Tref_);
    return props;
}

/* P [Pa], T [K] */
void hkf::dGMineral(double T, double P, double* G, double* S, double* a, double* b, double* c, double* V, int size,  double* dG)
{
    const double dT = T - Tref_;
    const double dP = P - Pref_;
    
    const double T_inv = 1. / T;
    const double dT2 = dT*dT;
    double A1 = (dT - T*log(T / Tref_));
    double A2 = -0.5*dT2;
    double A3 = -0.5*Tref2_inv_*T_inv*dT2;
    
    for (int i = 0; i < size; ++i)
    {
        dG[i] = G[i] - S[i] * dT + a[i] * A1 + b[i]*A2 + c[i]*A3 + V[i] * dP;
    }
}

/* NB: Assumes already KNOWN water properties at P [Pa] and T [K]. */
/* Note that we are actually using bar and as unit since Pressure is normalized to P_ref*/
void hkf::dGIons(double T, double P, double* G, double* S,
                 double* a1, double* a2, double* a3, double* a4,
                 double* c1, double* c2,
                 double* omega, double* Z, double* re_ref, int size,
                 int skip, double* dG, double* MV)
{
    const double dT = T - Tref_;
    
    const double Tp = T / Tref_;
    const double dP = P*Pref_inv_ - 1.0;
    
    const double C1 = dT - T*log(T / Tref_);
    const double C2 = -T / Theta_ / Theta_*log(Tp*(Tref_ - Theta_) / (T - Theta_)) - 1. / Theta_*(T - Tref_) / (Tref_ - Theta_);
    const double A1 = dP;
    const double A2 = log((Psi_ + P) / (Psi_ + Pref_));
    double A3 = 1.0 / (T - Theta_);
    const double A4 = A2*A3;
    const double Wref = Z_ref_ - EPS_->bornZ_ + Y_ref_*dT;
    const double W = -EPS_->bornZ_ - 1.0;
    
    // To calculate molar volume, the Born coefficient derivative is neglected
    //const double MV1 = 1.0 / (Psi_ + P); // @ah fix 27/2 2025
    const double MV1 = 1.0 / (Psi_*Pref_inv_ + P*Pref_inv_);
    const double MV2 = 1.0 / (T - Theta_);

    // const double Chat = 41.84 * 1e5 / UnitConversionFactors::cal2J_; //Pa*ml/J // @ah fix 27/2 2025

    const double Chat = 41.84*1e-3 / UnitConversionFactors::cal2J_; //bar*m^3/J
    
    A3 *= dP;
//    BORN_->born_f(T);
    BORN_->born_df(T, P);
    
    double ff;
    double Wi = 0.;
    double w_T, w_TT, w_P, Wi_2;

    for (int i = 0; i < size; ++i)
    {
        if (i != skip)
        {
            if (Z[i] == 0.)
            {
                ff = 0.;
                w_T = w_TT = w_P = 0.;
            }
            else
            {
//                BORN_->born(Z[i], re_ref[i], Wi);
                
                BORN_->born(Z[i], re_ref[i], Wi_2, w_T, w_TT, w_P);
                Wi = Wi_2;
 //               w_P = 0.;
                
                ff = (Wi-omega[i])*W;
            }
            dG[i] = G[i] - S[i] * dT + c1[i] * C1 + c2[i] * C2 + a1[i] * A1 + a2[i] * A2 + a3[i] * A3
            + a4[i] * A4 + ff + omega[i] * Wref;
            // Note that born_Q are in 1/Pa, while we switched back to bar to be consitent with original paper, same with w_P
            MV[i] = a1[i] + a2[i] * MV1 + (a3[i] + a4[i] * MV1)*MV2 - 1e5*omega[i]*EPS_->bornQ_-(EPS_->bornZ_+1)*w_P*1e5;
            MV[i] *= Chat;
            
            /* DEBUG*/
            /*
            double dg_delta= -S[i] * dT + c1[i] * C1 + c2[i] * C2 + a1[i] * A1 + a2[i] * A2 + a3[i] * A3
                + a4[i] * A4 + ff + omega[i] * Wref;
            std::cout << i << "eps_Tr_Pf" << 1 / Z_ref_ << " eps_T_P " << 1/EPS_->bornZ_ << std::endl;
            std::cout << i << " omega[i] = " << omega[i] << " J/mol " << omega[i] / 4.184 << " cal/mol "<< std::endl;
            std::cout << i << " re_ref[i] = " << re_ref[i] << std::endl;
            std::cout << i << " G = " << G[i] << " J/mol " << G[i] / 4.184 << " cal/mol" <<  G[i]/dG[i]*100 << " \%"<< std::endl;
            std::cout << i << " SdT = " << S[i] * dT << " J/mol " << S[i] * dT / 4.184 << " cal/mol " << S[i] * dT / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " c1[i] * C1 = " << c1[i] * C1 << " J/mol " << c1[i] * C1 / 4.184 << " cal/mol " << c1[i] * C1 / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " c2[i] * C2 = " << c2[i] * C2 << " J/mol " << c2[i] * C2 / 4.184 << " cal/mol " << c2[i] * C2 / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " a1[i] * A1 = " << a1[i] * A1 << " J/mol " << a1[i] * A1 / 4.184 << " cal/mol " << a1[i] * A1 / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " a2[i] * A2 = " << a2[i] * A2 << " J/mol " << a2[i] * A2 / 4.184 << " cal/mol " << a2[i] * A2 / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " a3[i] * A3 = " << a3[i] * A3 << " J/mol " << a3[i] * A3 / 4.184 << " cal/mol " << a3[i] * A3 / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " a4[i] * A4 = " << a4[i] * A4 << " J/mol " << a4[i] * A4 / 4.184 << " cal/mol " << a4[i] * A4 / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " (Wi-omega[i])*W = " << ff << " J/mol " << ff/ 4.184 << " cal/mol " << ff / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " Wi*W = " << Wi*W << " J/mol " << Wi*W / 4.184 << " cal/mol " << Wi*W / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << " omega[i] * Wref = " << omega[i] * Wref << " J/mol " << omega[i] * Wref / 4.184 << " cal/mol " << omega[i] * Wref / dG[i] * 100 << " \%" << std::endl;
            std::cout << i << "***** dG_delta = " << dg_delta << " J/mol " << dg_delta / 4.184 << " cal/mol " << std::endl;;

            std::cout << i << "***** dG = " << dG[i] << " J/mol " << dG[i] / 4.184 << " cal/mol " << std::endl;;
        */    
        }
    }
}

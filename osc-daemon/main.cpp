#include <iostream>
#include <fstream>
#include <string>

#include "osctypes.hpp"
#include "orbitalmechanics/axistransforms.hpp"


// get ASAT positions and velocities from following equations
// deltaM = (ASATflightTime - ReactionTime)*MeanOrbitalMotion(SatKOE)
// deltaE = MeanToTrue(deltaM, SatKOE)
// ASATp = [0, deltaE, h]
// ASATv = [V*cos(deltaE), V*sin(deltaE), 0]

// set up big fucking nasty array
// every combination between two orthogonal sets of deltaVs

// set up bigger fucking nastier array
// BFNA by 13 array

// for each combination
// perform burn of combination deltaV (Burn)
// find intercept values (InterceptCalcs)
// find if intercept is exo/endo atmospheric and viable
// some vile fucking logic here
// save any viable intercepts, ignore the rest

// find lowest deltaV burn for both exo and atmo

// do SMAD calculations for Casat
// plot graphs

// function [CaSatInterceptTime, CaSatAltAtTheta, ASAT_AltAtIntercept] = InterceptCalcs(KOE)
//     global deltaM mu EarthRadius ModelAccuracy ReactionTime h
//     %CASAT altitude equations
//     [dTheta,~]=MeanToTrue(deltaM, KOE(2));
//     if KOE(6)+dTheta>=360
//        dTheta=dTheta-360;
//     end
//     Function=@(v) 1./(1+KOE(2).*cosd(v)).^2;
//     if KOE(2)>0 && KOE(2)<1
//     [r, V] = keplerian2ijk(KOE(1), KOE(2), KOE(3), KOE(4), KOE(5), KOE(6)+dTheta, 'lonper', KOE(6)+dTheta);
//     elseif KOE(2)==0
//     [r, V] = keplerian2ijk(KOE(1), KOE(2), KOE(3), KOE(4), KOE(5), KOE(6)+dTheta, 'truelon', KOE(6)+dTheta);
//     else
//         CaSatInterceptTime=99999999;
//         CaSatAltAtTheta=999999999;
//         ASAT_AltAtIntercept=99999999;
//         return
//     end
//     AngularMomentum=norm(cross(r, V));
//     if KOE(6)>315 && KOE(6)+dTheta<45
//         dTheta=dTheta+360;
//     end
//     CaSatInterceptTime=abs(integral(Function, deg2rad(KOE(6)), deg2rad(KOE(6)+dTheta)))/(mu^2/AngularMomentum^3);
//     CaSatAltAtTheta=norm(r)-EarthRadius;
//     if ReactionTime/ModelAccuracy+(round(CaSatInterceptTime/ModelAccuracy))<=1000000
//         ASAT_AltAtIntercept=h(ReactionTime/ModelAccuracy+(round(CaSatInterceptTime/ModelAccuracy)));
//     else
//         ASAT_AltAtIntercept=10^9;
//     end

// end

// function KOE = Burn(DeltaVvnb, KOE) % burns by a[dVr dVt dVn] m / s impulse burn if KOE (2) == 0
//[r, t] = keplerian2ijk(KOE(1), KOE(2), KOE(3), KOE(4), KOE(5), KOE(6), 'truelon', KOE(6));
// else
//     [r, t] = keplerian2ijk(KOE(1), KOE(2), KOE(3), KOE(4), KOE(5), KOE(6), 'lonper', KOE(6));
// end
//     vnb2ijk = [t./ norm(t) cross(r, t)./ norm(cross(r, t)) cross(t, cross(r, t))./ norm(cross(t, cross(r, t)))];
// DeltaVijk = vnb2ijk * DeltaVvnb'; t = t + DeltaVijk;
//                           [KOE(1), KOE(2), KOE(3), KOE(4), KOE(5), KOE(6), trulon] = ijk2keplerian(r, t);
// if
//     ~isfinite(KOE(6)) == 1 KOE(6) = trulon;
// end if ~isfinite(KOE(4)) == 1 || ~isfinite(KOE(5)) == 1 KOE(4) = 0;f
// KOE(5) = 0;
// end
//     end

//         function n = MeanOrbitalMotion(KOE)
//             global mu
//                 n = rad2deg(sqrt(mu / pow3(KOE(1))));
// end
        // needs to create a KOE object with these parameters
    //  Sat_Alt=250000; Sat_Ecc=0;      Sat_Inc=0;
    //  Sat_RAAN=0;     Sat_ArgPer=0;   Sat_MeaAno=0;
    //  Sat_SMA=EarthRadius+Sat_Alt;
    //  SatKOE=[Sat_SMA Sat_Ecc Sat_Inc Sat_RAAN Sat_ArgPer Sat_MeaAno];

    // constexpr const double sat_ecc = 0;
    // constexpr const double sat_inc = 0;
    // constexpr const double sat_RAAN = 0;
    // constexpr const double sat_ArgPer = 0;
    // constexpr const double sat_MeanAno = 0;
    // constexpr const double sat_alt = 250000;
    // constexpr const double semi_maj_ax = sat_alt + EarthRadius;

    //can i get away with doing this?

osc::orbParam Burn(double v, double n, double b, osc::orbParam KOE){
    using namespace osc;
    vnb dV;
    dV.vVNB = {v, n, b};

    pcs PCSposvel = KOEtoPCS(KOE);
    eci ECIposvel = PCStoECI(KOE, PCSposvel);
    eci iECI = VNBtoECI(ECIposvel, dV);
    orbParam iKOE = ECItoKOE(iECI);
    return iKOE;
}

osc::vec3 InterceptCalcs(double deltaM, osc::orbParam KOE){
    using namespace osc;
    double dTheta = KOE.meanToTrue(deltaM);
    if (KOE.truAnom+dTheta>=2 * M_PI){dTheta=dTheta-2 * M_PI;}

    pcs PCSposvel = KOEtoPCS(KOE);
    eci ECIposvel = PCStoECI(KOE, PCSposvel);

    //i cannot remember why i added this
    if (KOE.truAnom>7*M_PI_4 && KOE.truAnom+dTheta<M_PI_4){dTheta=dTheta+2*M_PI;}

    vec3 h = ECIposvel.rIJK.cross(ECIposvel.vIJK);//angular momentum
    double h_abs = h.mag();

    double CaSatInterceptTime = KOE.RK4(dTheta, 0.001)/(pow2(mu)/pow3(h_abs));
    double CaSatAltAtTheta = ECIposvel.rIJK.mag()-EarthRadius;
    double ASAT_AltAtIntercept; //needs to read that .csv file from before
    return {CaSatInterceptTime, CaSatAltAtTheta, ASAT_AltAtIntercept};
}

int main(int, char **)
{



    using namespace std;
    ofstream outputFile;
    outputFile.open("AsatM_Out.csv");

    using namespace osc;
    int LRCT;
    int HRCT;
    double ASATAltAtLRCT;
    struct orbParam SatKOE = {250000+EarthRadius, 0, 0, 0, 0, 0};
    struct orbParam CaSatKOE = {250000+EarthRadius, 0, 0, 0, 0, 0};
    double deltaM = (t_f-ReactionTime)*SatKOE.MeanOrbitalMotion();
    double deltaE = SatKOE.meanToTrue(SatKOE.ecc);

    // defender satellite initial keplerian orbital elements
    //  CaSat_Alt_i=250000; CaSat_Ecc_i=0;      CaSat_Inc_i=0;
    //  CaSat_RAAN_i=0;     CaSat_ArgPer_i=0;   CaSat_MeaAno_i=0;
    //  CaSat_SMA_i=EarthRadius+CaSat_Alt_i;
    //  CaSatKOEi=[CaSat_SMA_i CaSat_Ecc_i CaSat_Inc_i CaSat_RAAN_i CaSat_ArgPer_i CaSat_MeaAno_i];

    // ASAT-M Model Values
    double mo = ((Mo / std::exp(((V_f / g) + t_f) / 280)) - Mo) / (-t_f);      // fuel mass flow Rate
    double p_f = 1 - ((Mo - mo * t_f) / Mo) * (log(Mo / (Mo - mo * t_f)) + 1); // range function p(f)
    double TWR = (g * pow2(Isp) * p_f) / (h_f + 0.5 * g * pow2(t_f));          // thrust to weight Ratio

    // ASAT-M Modelling with time
    // t=ModelAccuracy:ModelAccuracy:t_f;              //sets up matrix for all times
   for (double t = 0; t < t_f; t+=ModelAccuracy){

    double p = 1 - ((Mo - mo * t) / Mo) * (log(Mo / (Mo - mo * t)) + 1); // range function p(f) again
    double h = ((g * pow2(Isp) / TWR) * p - (0.5 * g * pow2(t)));        // altitude as function of time
    double V = g * (Isp * log(Mo / (Mo - mo * t)) + 1);                  // velocity as function of time
    double h_c = h + (pow2(V)) / (2 * g);                                // culmination alt as func of time
    if (h_c>100000){
        LRCT=t;
        ASATAltAtLRCT=h;};
    if (h>(SatKOE.sma-EarthRadius)){HRCT=t;};

    lla LLApos = {0, deltaE, h};
    ecef asatECEF = LLAtoECEF(LLApos);
    asatECEF.vXYZ = {V*cos(deltaE), V*sin(deltaE), 0};

    outputFile << h << "," << V << "," <<  h_c << "," << asatECEF.rXYZ[0] << "," << asatECEF.rXYZ[1] << "," << asatECEF.rXYZ[2] << "," <<  asatECEF.vXYZ[0] << "," << asatECEF.vXYZ[1] << "," << asatECEF.vXYZ[2] << std::endl;

   }
   outputFile.close();

    outputFile.open("CaSatM_Out.csv");
    for (double dVv = -150; dVv<=-50; dVv+=1){
        for (double dVb = -1250; dVb<=-250; dVb+=1){
            orbParam iKOE = Burn(dVv, 0, dVb, SatKOE);
            vec3 InterceptOut = InterceptCalcs(deltaM, SatKOE);
            outputFile << dVv << "," << dVb << "," <<  sqrt(pow2(dVv)+pow2(dVb)) << "," << iKOE.sma << "," << iKOE.ecc << "," << iKOE.inc << "," <<  iKOE.aop << "," << iKOE.asc << "," << iKOE.truAnom << "," << InterceptOut[0] << InterceptOut[1]/1000 << InterceptOut[2]/1000 << std::endl;
        }
        std::cout << dVv; //debug
    };
}
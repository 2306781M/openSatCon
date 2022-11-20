#include <fstream>
#include <iostream>
#include <string>

#include "orbitalmechanics/axistransforms.hpp"
#include "osctypes.hpp"

osc::orbParam Burn(double v, double n, double b, osc::orbParam KOE) {
  using namespace osc;
  vnb dV;
  dV.vVNB = {v, n, b};

  // this doesnt need to be done every time, nothing actually changes until
  // ECIdV
  pcs PCSposvel = KOEtoPCS(KOE);
  eci ECIposvel = PCStoECI(KOE, PCSposvel);
  eci ECIdV = VNBtoECI(ECIposvel, dV);
  // std::cout << dV.vVNB[0] << ",    "<< dV.vVNB[1] << ",    "<< dV.vVNB[2] <<
  // ",    "<<std::endl;
  eci iECI;
  iECI.rIJK = ECIposvel.rIJK;
  iECI.vIJK = ECIposvel.vIJK.operator+(ECIdV.vIJK);
  orbParam iKOE = ECItoKOE(iECI);
  return iKOE;
}

osc::vec3 InterceptCalcs(double deltaM, osc::orbParam KOE) {
  using namespace osc;
  double dTheta = KOE.meanToTrue(deltaM);
  KOE.truAnom+=dTheta;
  pcs PCSposvel = KOEtoPCS(KOE);
  // std::cout << PCSposvel.rPCS[0] << ",   " << PCSposvel.rPCS[1] << ",   " <<
  // PCSposvel.rPCS[2] << std::endl;
  eci ECIposvel = PCStoECI(KOE, PCSposvel);
  KOE.truAnom-=dTheta;
  double CaSatInterceptTime = KOE.CompositeTrapezoid(dTheta);
  double CaSatAltAtTheta = ECIposvel.rIJK.mag() - EarthRadius;
  // std::cout << ECIposvel.rIJK[0] << ",    " << ECIposvel.rIJK[1] << ",    "
  //           << ECIposvel.rIJK[2] << ",    " << ECIposvel.rIJK.mag()
  //           << std::endl;
  double ASAT_AltAtIntercept; // needs to read that .csv file from before
  return {CaSatInterceptTime, CaSatAltAtTheta, ASAT_AltAtIntercept};
}

int main(int, char **) {
  using namespace std;
  ofstream outputFile;
  outputFile.open("AsatM_Out.csv", std::ofstream::trunc);

  using namespace osc;
  double h;
  double V;
  double LRCT;
  double HRCT;
  double deltaM;
  double ASATAltAtLRCT;
  struct orbParam SatKOE = {250000 + EarthRadius, 0, 0, 0, 0, 0};
  struct orbParam CaSatKOE = {250000 + EarthRadius, 0, 0, 0, 0, 0};

  // defender satellite initial keplerian orbital elements
  //  CaSat_Alt_i=250000; CaSat_Ecc_i=0;      CaSat_Inc_i=0;
  //  CaSat_RAAN_i=0;     CaSat_ArgPer_i=0;   CaSat_MeaAno_i=0;
  //  CaSat_SMA_i=EarthRadius+CaSat_Alt_i;
  //  CaSatKOEi=[CaSat_SMA_i CaSat_Ecc_i CaSat_Inc_i CaSat_RAAN_i CaSat_ArgPer_i
  //  CaSat_MeaAno_i];

  // ASAT-M Model Values
  double mo = ((Mo / std::exp(((V_f / g) + t_f) / 280)) - Mo) /
              (-t_f); // fuel mass flow Rate
  double p_f = 1 - ((Mo - mo * t_f) / Mo) *
                       (log(Mo / (Mo - mo * t_f)) + 1); // range function p(f)
  double TWR = (g * pow2(Isp) * p_f) /
               (h_f + 0.5 * g * pow2(t_f)); // thrust to weight Ratio

  // ASAT-M Modelling with time
  // t=ModelAccuracy:ModelAccuracy:t_f;              //sets up matrix for all
  // times
  for (double t = 0; t < t_f; t += ModelAccuracy) {

    double p = 1 - ((Mo - mo * t) / Mo) * (log(Mo / (Mo - mo * t)) +
                                           1); // range function p(f) again
    double h = ((g * pow2(Isp) / TWR) * p -
                (0.5 * g * pow2(t))); // altitude as function of time
    double V =
        g * (Isp * log(Mo / (Mo - mo * t)) + 1); // velocity as function of time
    double h_c = h + (pow2(V)) / (2 * g); // culmination alt as func of time
    if (h_c > 100000) {
      LRCT = t;
      ASATAltAtLRCT = h;
    };
    if (h <= (SatKOE.sma - EarthRadius)) {
      HRCT = t;
    };

    outputFile << h << "," << V << "," << h_c << std::endl;
  }
  outputFile.close();

  outputFile.open("CaSatMposVel_Out.csv", std::ofstream::trunc);
  deltaM = (HRCT - ReactionTime) * SatKOE.MeanOrbitalMotion();
  double deltaE = SatKOE.meanToTrue(deltaM);
  // for (int i=0; i<0; i++){
  //     lla LLApos = {0, deltaE, h};
  //     ecef asatECEF = LLAtoECEF(LLApos);
  //     //std::cout << LLApos.lat << "  ,   " << LLApos.lon << "  ,   " <<
  //     LLApos.alt << std::endl;
  //     //std::cout << asatECEF.rXYZ[0] << "  ,   " << asatECEF.rXYZ[1] << "  ,
  //     " << asatECEF.rXYZ[2] << std::endl << std::endl; asatECEF.vXYZ =
  //     {V*cos(deltaE), V*sin(deltaE), 0}; outputFile << asatECEF.rXYZ[0] <<
  //     "," << asatECEF.rXYZ[1] << "," << asatECEF.rXYZ[2] << "," <<
  //     asatECEF.vXYZ[0] << "," << asatECEF.vXYZ[1] << "," << asatECEF.vXYZ[2]
  //     << std::endl;
  // }
  outputFile.close();

  outputFile.open("CaSatM_Out.csv", std::ofstream::trunc);
  for (double dVb = -350; dVb <= -50; dVb += 10) {
    for (double dVv = -1250; dVv <= -250; dVv += 10) {
      orbParam iKOE = Burn(dVv, 0, dVb, SatKOE);
      vec3 InterceptOut = InterceptCalcs(deltaM, iKOE);
      outputFile << dVv << "," << dVb << "," << sqrt(pow2(dVv) + pow2(dVb))
                 << "," << iKOE.sma << "," << iKOE.ecc << "," << iKOE.inc << ","
                 << iKOE.aop << "," << iKOE.asc << "," << iKOE.truAnom << ","
                 << InterceptOut[0] << "," << InterceptOut[1] / 1000 << ","
                 << InterceptOut[2] / 1000 << std::endl;
    }
  };
  std::cout << std::endl << "complete" << std::endl;
}
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "orbitalmechanics/axistransforms.hpp"
#include "osctypes.hpp"
#include <gcem.hpp>
#include <oneapi/tbb.h>

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

osc::vec3 InterceptCalcs(double deltaM, osc::orbParam KOE, double AAALRCT, double HRCT,
                         double SatAlt, double dVv, double dVb, double h[]) {
  using namespace osc;
  double dTheta = KOE.meanToTrue(deltaM);
  double CaSatInterceptTime;
  double ASAT_AltAtIntercept;
  double CaSatAltAtTheta;
  KOE.truAnom += dTheta;
  pcs PCSposvel = KOEtoPCS(KOE);
  // PCSposvel.rPCS[2] << std::endl;
  eci ECIposvel = PCStoECI(KOE, PCSposvel);
  KOE.truAnom -= dTheta;
  CaSatInterceptTime = KOE.CompositeTrapezoid(dTheta);
  if (CaSatInterceptTime+ReactionTime<HRCT){
    double CaSatAltAtTheta = ECIposvel.rIJK.mag() - EarthRadius;
  ASAT_AltAtIntercept = h[int(round(((ReactionTime + CaSatInterceptTime) / ModelAccuracy)))]; 
  if (abs(CaSatAltAtTheta - ASAT_AltAtIntercept) < 10000) {
    if (CaSatAltAtTheta < AAALRCT + 5000) {
      //std::cout << "Successful Endoatmospheric Intercept at " << dVv
                //<< "m/s prograde, " << dVb << "m/s radial.\n"
                //<< "Intercept Altitude: " << CaSatAltAtTheta << "m\n\n";
    } else if (ASAT_AltAtIntercept > 100000 &&
               ASAT_AltAtIntercept < (SatAlt - 50000) &&
               abs(CaSatAltAtTheta - ASAT_AltAtIntercept) < 1000) {
      //std::cout << "Successful Exoatmospheric Intercept at " << dVv
                //<< "m/s prograde, " << dVb << "m/s radial.\n"
                //<< "Intercept Altitude: " << CaSatAltAtTheta << "m\n\n";
    }
  }
  }
  else {CaSatAltAtTheta = 9e10; ASAT_AltAtIntercept = 9e10;}
  return {CaSatInterceptTime, CaSatAltAtTheta, ASAT_AltAtIntercept};
}

int main() {
  using namespace std;
  ofstream outputFile;
  outputFile.open("AsatM_Out.csv", std::ofstream::trunc);

  using namespace osc;

  auto h = std::make_unique<std::array<double, 168000000>>();
  double V=0;
  double LRCT;
  double HRCT;
  double deltaM;
  double ASATAltAtLRCT=0;
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
  int i = 0;
  for (double t = 0; t < t_f; t += ModelAccuracy) {

    double p = 1 - ((Mo - mo * t) / Mo) * (log(Mo / (Mo - mo * t)) +
                                           1); // range function p(f) again
    (*h)[i] = ((g * pow2(Isp) / TWR) * p -
               (0.5 * g * pow2(t))); // altitude as function of time
    V =
        g * (Isp * log(Mo / (Mo - mo * t)) + 1); // velocity as function of time
    double h_c =
        (*h)[i] + (pow2(V)) / (2 * g); // culmination alt as func of time
    if (h_c < 100000) {
      LRCT = t;
      ASATAltAtLRCT = (*h)[i];
    };
    if ((*h)[i] <= (SatKOE.sma - EarthRadius)) {
      HRCT = t;
    };

    outputFile << (*h)[i]  << "," << V << "," << h_c << std::endl;
    i++;
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

  /**
   * Simple Helper to linspace a vector
   */
  auto linspace = [](auto &vec, auto start, auto stop, size_t n_values) {
    vec.reserve(n_values);
    vec.resize(n_values);
    using vec_value_t = std::remove_cvref_t<decltype(vec)>::value_type;
    vec_value_t increment = static_cast<vec_value_t>((stop - start) / n_values);
    std::generate(vec.begin(), vec.end(),
                  [n = 0, &increment, &start]() mutable -> vec_value_t {
                    return (n++ * increment) + start;
                  });
  };

  outputFile.open("CaSatM_Out.csv", std::ofstream::trunc);
  outputFile.close();

  auto do_calculation = [&](double dVb, double dVv) {
    orbParam iKOE = Burn(dVv, 0, dVb, SatKOE);
    vec3 InterceptOut =
        InterceptCalcs(deltaM, iKOE, ASATAltAtLRCT, HRCT, SatKOE.sma - EarthRadius,
                       dVv, dVb, h.get()->data());
    // outputFile << dVv << "," << dVb << "," << sqrt(pow2(dVv) + pow2(dVb)) <<
    // ","
    //            << iKOE.sma << "," << iKOE.ecc << "," << iKOE.inc << ","
    //            << iKOE.aop << "," << iKOE.asc << "," << iKOE.truAnom << ","
    //            << InterceptOut[0] << "," << InterceptOut[1] / 1000 << ","
    //            << InterceptOut[2] / 1000 << std::endl;
  };

  std::vector<double> dVbs{};
  std::vector<double> dVvs{};
  // Adjust linspece values based on needs
  linspace(dVbs, -5000.0, 5000.0, 1000);
  linspace(dVvs, -5000.0, 5000.0, 1000);
  tbb::parallel_for(
      tbb::blocked_range2d<double>(0, dVbs.size(), 0, dVbs.size()),
      [&](const tbb::blocked_range2d<double> &range) {
        for (int i = range.rows().begin(), i_end = range.rows().end();
             i < i_end; ++i) {
          for (int j = range.cols().begin(), j_end = range.cols().end();
               j < j_end; ++j) {
            do_calculation(dVbs[i], dVvs[i]);
          }
        }
      });
  std::cout << "\n\n" << "COMPLETE" << "\n\n";
}
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "orbitalmechanics/axistransforms.hpp"
#include "osctypes.hpp"
#include <gcem.hpp>
#include <mutex>
#include <oneapi/tbb.h>
#include <optional>
#include <type_traits>

osc::eci LambertsProblem(osc::vec3 r1, osc::vec3 r2, double t) {
  using namespace osc;

  vec3 c = r2 - r1;

  vec3 t1Unit, t2Unit;
  double r1Mag = r1.mag();
  double r2Mag = r2.mag();
  double cMag = c.mag();

  double s = 0.5 * (r1Mag + r2Mag + cMag);
  double x;
  vec3 r1Unit = r1.unit();

  vec3 r2Unit = r2.unit();

  vec3 hUnit = r1Unit.cross(r2Unit).unit();
  double lambda = sqrt(1.0 - cMag / s);

  if ((hUnit[2]) < 0.0) {
    lambda = -lambda;
    t1Unit = r1Unit.cross(hUnit);
    t2Unit = r2Unit.cross(r2Unit);
  } else {
    t1Unit = hUnit.cross(r1Unit);
    t2Unit = hUnit.cross(r2Unit);
  }
  double lambda2 = pow2(lambda);
  double lambda3 = pow3(lambda);

  double T = sqrt(2.0 * mu / pow3(s)) * t;

  double T0 = acos(lambda) + lambda * sqrt(1 - lambda2);

  double T1 = (2.0 / 3.0) * (1.0 - lambda3);

  if (T >= T0) {
    x = pow(T0 / T, 2.0 / 3.0) - 1.0;
  } else if (T < T1) {
    x = 2.5 * ((T1 * (T1 - T)) / (T * (1.0 - lambda2 * lambda3))) + 1.0;
  } else {
    x = pow((T / T0), 0.69314718055994529 / log(T1 / T0)) - 1.0;
  }
  int i = 0;
  int iMax = 15; // idk what to put
  double err = 1.0;
  double eps = 1e-5;
  double xn = 0.0;
  double tof = 0.0, delta = 0.0, deriv1 = 0.0, deriv2 = 0.0, deriv3 = 0.0;
  while ((err > eps) && (i < iMax)) {

    double battin = 0.01;
    double lagrange = 0.2;
    double distance = abs(x - 1.0);
    if (distance < lagrange && distance > battin) {
      double a = 1.0 / (1.0 - pow2(x));
      if (a > 0.0) {
        double alpha = 2.0 * acos(x);
        double beta = 2.0 * asin(sqrt(lambda2 / a));
        if (lambda < 0.0) {
          beta = -beta;
        }
        tof =
            ((a * sqrt(a) * ((alpha - sin(alpha)) - (beta - sin(beta)))) / 2.0);
      } else {
        double alpha = 2.0 * acosh(x);
        double beta = 2.0 * asinh(sqrt(-lambda * lambda / a));
        if (lambda < 0.0) {
          beta = -beta;
        }
        tof = ((-a * sqrt(-a) * ((beta - sinh(beta)) - (alpha - sinh(alpha)))) /
               2.0);
      }
    }
    double E = pow2(x) - 1.0;
    double rho = abs(E);
    double z = sqrt(1.0 + lambda2 * E);
    if (distance < battin) {
      double eta = z - lambda * x;
      double S1 = 0.5 * (1.0 - lambda - x * eta);
      double Sj = 1.0;
      double Cj = 1.0;
      double err = 1.0;
      double Cj1 = 0.0;
      double Sj1 = 0.0;
      int j = 0;
      while (err > 1e-11) {
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * S1 / (j + 1.0);
        Sj1 = Sj + Cj1;
        err = abs(Cj1);
        Sj = Sj1;
        Cj = Cj1;
        j++;
      }
      double Q = 4.0 / 3.0 * Sj;
      tof = (pow3(eta) * Q + 4.0 * lambda * eta) / 2.0;
    } else {
      double y = sqrt(rho);
      double g = x * z - lambda * E;
      double d = 0.0;
      if (E < 0.0) {
        d = acos(g);
      } else {
        double f = y * (z - lambda * x);
        d = log(f + g);
      }
      tof = (x - lambda * z - d / y) / E;
    }

    double umx2 = 1.0 - pow2(x);
    double invumx2 = 1.0 / umx2;
    double y = sqrt(1 - lambda2 * umx2);
    double y2 = pow2(y);
    double y3 = pow3(y);
    double deriv1 = invumx2 * (3.0 * T * x - 2.0 + 2.0 * lambda3 * x / y);
    double deriv2 = invumx2 * (3.0 * T + 5.0 * x * deriv1 +
                               2.0 * (1.0 - lambda2) * lambda3 / y3);
    double deriv3 =
        invumx2 * (7.0 * x * deriv2 + 8.0 * deriv1 -
                   6.0 * (1.0 - lambda2) * lambda2 * lambda3 * x / y3 / y2);
    delta = tof - T;
    double deriv12 = pow2(deriv1);
    xn = x - delta * (deriv12 - delta * deriv2 / 2.0) /
                 (deriv1 * (deriv12 - delta * deriv2) +
                  deriv3 * delta * delta / 6.0);
    err = abs(x - xn);
    x = xn;
    i++;
  }

  double gamma = sqrt(0.5 * mu * s);
  double rho = (r1Mag - r2Mag) / cMag;
  double sigma = sqrt(1.0 - pow2(rho));

  double y = sqrt(1.0 - lambda2 + lambda2 * pow2(x));
  double Vr1 = gamma * ((lambda * y - x) - rho * (lambda * y + 1.0)) / r1Mag;
  double Vr2 = -gamma * ((lambda * y - x) + rho * (lambda * y + 1.0)) / r2Mag;
  double Vt1 = gamma * sigma * (y + lambda * x) / r1Mag;
  double Vt2 = gamma * sigma * (y + lambda * x) / r2Mag;
  vec3 v1 = r1Unit * Vr1 + t1Unit * Vt1;
  vec3 v2 = r2Unit * Vr2 + t2Unit * Vt2;

  eci eciOut{r1, v1};

  return (eciOut);
}

osc::orbParam Burn(double v, double n, double b, osc::orbParam KOE) {
  using namespace osc;
  vnb dV;
  dV.vVNB = {v, n, b};

  // this doesnt need to be done every time, nothing actually changes until
  // ECIdV
  pcs PCSposvel = KOEtoPCS(KOE);
  eci ECIposvel = PCStoECI(KOE, PCSposvel);
  eci ECIdV = VNBtoECI(ECIposvel, dV);

  eci iECI;
  iECI.rIJK = ECIposvel.rIJK;
  iECI.vIJK = ECIposvel.vIJK.operator+(ECIdV.vIJK);
  orbParam iKOE = ECItoKOE(iECI);
  return iKOE;
}

std::optional<osc::InterceptOut>
InterceptCalcs(double deltaM, osc::orbParam KOE, double AAALRCT, double HRCT,
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
  if (CaSatInterceptTime + ReactionTime < HRCT) {
    double CaSatAltAtTheta = ECIposvel.rIJK.mag() - EarthRadius;
    ASAT_AltAtIntercept =
        h[int(round(((ReactionTime + CaSatInterceptTime) / ModelAccuracy)))];
    if (abs(CaSatAltAtTheta - ASAT_AltAtIntercept) < 10000) {
      if (CaSatAltAtTheta < AAALRCT + 5000) {
        InterceptOut EndoAtmoIntercept{
            dVv,
            dVb,
            sqrt(pow2(dVv) + pow2(dVb)),
            KOE,
            {CaSatInterceptTime, CaSatAltAtTheta, ASAT_AltAtIntercept}};
        return EndoAtmoIntercept;
      } else if (ASAT_AltAtIntercept > 100000 &&
                 ASAT_AltAtIntercept < (SatAlt - 50000) &&
                 abs(CaSatAltAtTheta - ASAT_AltAtIntercept) < 1000) {
        InterceptOut ExoAtmoIntercept{
            dVv,
            dVb,
            sqrt(pow2(dVv) + pow2(dVb)),
            KOE,
            {CaSatInterceptTime, CaSatAltAtTheta, ASAT_AltAtIntercept}};
        return ExoAtmoIntercept;
        // std::cout << "Successful Exoatmospheric Intercept at " << dVv
        //<< "m/s prograde, " << dVb << "m/s radial.\n"
        //<< "Intercept Altitude: " << CaSatAltAtTheta << "m\n\n";
      }
    }
  }
  // else {CaSatAltAtTheta = 9e10; ASAT_AltAtIntercept = 9e10;}
  // return {CaSatInterceptTime, CaSatAltAtTheta, ASAT_AltAtIntercept};
  return std::nullopt;
}

int main() {
  using namespace std;
  ofstream outputFile;
  outputFile.open("AsatM_Out.csv", std::ofstream::trunc);

  using namespace osc;

  auto h = std::make_unique<std::array<double, 168000>>();
  auto V = std::make_unique<std::array<double, 168000>>();
  auto r2 = std::make_unique<std::array<osc::vec3, 168000>>();
  double LRCT;
  double HRCT;
  double deltaM;
  double ASATAltAtLRCT = 0;
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
    (*V)[i] =
        g * (Isp * log(Mo / (Mo - mo * t)) + 1); // velocity as function of time
    double h_c =
        (*h)[i] + (pow2((*V)[i])) / (2 * g); // culmination alt as func of time
    if (h_c < 100000) {
      LRCT = t;
      ASATAltAtLRCT = (*h)[i];
    };
    if ((*h)[i] <= (SatKOE.sma - EarthRadius)) {
      HRCT = t;
    };

    outputFile << (*h)[i] << "," << (*V)[i] << "," << h_c << "\n";
    i++;
  }
  outputFile.close();

  outputFile.open("CaSatMposVel_Out.csv", std::ofstream::trunc);
  deltaM = (HRCT - ReactionTime) * SatKOE.MeanOrbitalMotion();
  double deltaE = SatKOE.meanToTrue(deltaM);
  for (int i = 0; i < 168000; i++) {
    lla LLApos = {0, deltaE, (*h)[i]};
    ecef asatECEF = LLAtoECEF(LLApos);
    (*r2)[i] = asatECEF.rXYZ;

    asatECEF.vXYZ = {(*V)[i] * cos(deltaE), (*V)[i] * sin(deltaE), 0};
    outputFile << asatECEF.rXYZ[0] << "," << asatECEF.rXYZ[1] << ","
               << asatECEF.rXYZ[2] << "," << asatECEF.vXYZ[0] << ","
               << asatECEF.vXYZ[1] << "," << asatECEF.vXYZ[2] << "\n";
  }
  outputFile.close();

  /*alternative main
  pcs PCSposvel = KOEtoPCS(SatKOE);
  eci ECIposvel = PCStoECI(SatKOE, PCSposvel);

  */
  pcs PCSposvel = KOEtoPCS(SatKOE);
  eci ECIposvel = PCStoECI(SatKOE, PCSposvel);

  outputFile.open("LambertSolutions.csv", std::ofstream::trunc);
  for (int i = ReactionTime / ModelAccuracy; i < HRCT / ModelAccuracy; i++) {
    double InterceptTime = i * ModelAccuracy - ReactionTime;
    eci outECI = LambertsProblem(ECIposvel.rIJK, (*r2)[i], InterceptTime);
    vec3 deltaV = outECI.vIJK - ECIposvel.vIJK;
    outputFile << (*h)[i] << "," << InterceptTime << "," << deltaV[0] << ","
               << deltaV[1] << "," << deltaV[2] << "," << deltaV.mag() << "\n";
  }
  outputFile.close();
  std::cout << "\n"
            << "LAMBERT SOLUTIONS COMPLETE"
            << "\n";
  /**
   * Simple Helper to linspace a vector
   */
  std::mutex output_mut;
  std::vector<osc::InterceptOut> output;
  auto linspace = [](auto &vec, auto start, auto stop, size_t n_values) {
    vec.reserve(n_values);
    vec.resize(n_values);
    using vec_value_t =
        std::remove_reference_t<std::remove_cv_t<decltype(vec)>>::value_type;
    vec_value_t increment = static_cast<vec_value_t>((stop - start) / n_values);
    std::generate(vec.begin(), vec.end(),
                  [n = 0, &increment, &start]() mutable -> vec_value_t {
                    return (n++ * increment) + start;
                  });
  };

  outputFile.open("CaSatM_Out.csv", std::ofstream::trunc);
  outputFile.close();

  auto DoCalculation = [&](double dVb, double dVv) -> void {
    orbParam iKOE = Burn(dVv, 0, dVb, SatKOE);
    auto InterceptOut =
        InterceptCalcs(deltaM, iKOE, ASATAltAtLRCT, HRCT,
                       SatKOE.sma - EarthRadius, dVv, dVb, h.get()->data());
    if (!InterceptOut.has_value()) {
      return;
    }
    {
      std::lock_guard<std::mutex> guard(output_mut);
      output.push_back(InterceptOut.value());
    }
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
  linspace(dVbs, -5000.0, 5000.0, 100);
  linspace(dVvs, -5000.0, 5000.0, 100);
  tbb::parallel_for(
      tbb::blocked_range2d<double>(0, dVbs.size(), 0, dVbs.size()),
      [&](const tbb::blocked_range2d<double> &range) {
        for (int i = range.rows().begin(), i_end = range.rows().end();
             i < i_end; ++i) {
          for (int j = range.cols().begin(), j_end = range.cols().end();
               j < j_end; ++j) {
            DoCalculation(dVbs[i], dVvs[i]);
          }
        }
      });
  std::cout << "\n\n"
            << "COMPLETE"
            << "\n\n";
}
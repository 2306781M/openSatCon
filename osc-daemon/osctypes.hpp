#ifndef OSCTYPES_H
#define OSCTYPES_H

#include <array>
#include <cassert>
#include <math.h>
#include <vector>

namespace osc {
constexpr const double g = 9.80665;
constexpr const int h_f = 283000;                 // final height of ASAT-M
constexpr const int t_f = 168;                    // time to reach final height
constexpr const int V_f = 3400;                   // final velocity of ASAT-M
constexpr const int Mo = 19000;                   // initial mass of ASAT-M
constexpr const double Lambda = 19 / (19 - 16.7); // Payload fraction values
constexpr const int Isp = 280; // standard value for solid fuel motor
constexpr const double ModelAccuracy = 0.001; // modelling value
constexpr const double mu = 3.986012e14;      // keeping large values to avoid
constexpr const int EarthRadius = 6371000;    // unit changing mid-model
constexpr const int AltAtmosphere = 84852 + EarthRadius;
constexpr const int ReactionTime = 20;

// Inline helper functions
/// Inline function to square an argument
constexpr double pow2(double arg) { return std::pow(arg, 2); }

/// Inline function to cube an argument
constexpr double pow3(double arg) { return std::pow(arg, 3); }

constexpr double deg2rad(double arg) { return arg * M_PI / 180; }
constexpr double rad2deg(double arg) { return arg / M_PI * 180; }

// Common data types
/** \struct vec3
\brief vector3
Vector3 datatype
*/
struct vec3 {
  /// Class variable initialiser
  std::array<double, 3> data = {0, 0, 0};

  // Initialisers
  constexpr vec3(std::array<double, 3> initArray) : data(initArray) {}

  constexpr vec3(double i, double j, double k) : data({i, j, k}) {}

  constexpr vec3() {}

  // Operator overloads
  /// Addition Compatibility between two vec3s
  constexpr vec3 operator+(const vec3 &rhs) const noexcept {
    return vec3(data[0] + rhs[0], data[1] + rhs[1], data[2] + rhs[2]);
  }

  // vec3 operatorplus2 (vec3 term1, vec3 term2) { return vec3(data[0] +
  // term1[0] + term2[0], data[1] + term1[1] + term2[1], data[2] + term1[2] +
  // term2[2]); }
  /// Subtraction Compatibility between two vec3s
  constexpr vec3 operator-(const vec3 &rhs) const noexcept {
    return vec3(data[0] - rhs[0], data[1] - rhs[1], data[2] - rhs[2]);
  }
  /// Multiplication Compatibility between two vec3s
  constexpr vec3 operator*(double rhs) const noexcept {
    return vec3(data[0] * rhs, data[1] * rhs, data[2] * rhs);
  }
  /// Division Compatibility between two vec3s
  constexpr vec3 operator/(double rhs) const noexcept {
    return vec3(data[0] / rhs, data[1] / rhs, data[2] / rhs);
  }
  /// Adding usage of [] operator
  constexpr double operator[](size_t i) const noexcept {
    assert(i >= 0 && i < 3 && "invalid Index");
    return data[i];
  }

  /// Implicit type conversion
  constexpr operator std::array<double, 3>() const noexcept { return data; };

  // Member functions
  /** \fn dot(arg)
  @param[in] arg
  @param[out] vec3
  does the dot product with the vec3 and another vec3
  */
  constexpr double dot(const vec3 &arg) const noexcept {
    /*
    Performs the vector dot product with the argument
    */
    return (data[0] * arg[0] + data[1] * arg[1] + data[2] * arg[2]);
  }

  /** \fn cross(arg)
  @param[in] arg
  @param[out] vec3
  does the cross product with the vec3 and another vec3
  */
  constexpr vec3 cross(const vec3 &arg) const noexcept {
    /*
    Performs the vector cross product with the argument
    */
    return vec3(data[1] * arg[2] - data[2] * arg[1],
                data[2] * arg[0] - data[0] * arg[2],
                data[0] * arg[1] - data[1] * arg[0]);
  }

  /** \fn mag()
  @param[out] double
  calculates the magnitude of the vec3
  */
  constexpr double mag() const noexcept { return sqrt(dot(*this)); }
  /** \fn unit()
  @param[out] vec3
  calculates the unit vector of the vec3
  */
  constexpr vec3 unit() const noexcept { return vec3(data) / mag(); };
};

/** \struct position
\brief x,y,z position
position struct to define an object's position
*/
struct position {
  // Class variable initialisers
  double x = 0;
  double y = 0;
  double z = 0;

  // Initialisers
  constexpr position(double initX, double initY, double initZ) noexcept
      : x(initX), y(initY), z(initZ) {}
  constexpr position(std::array<double, 3> initPos) noexcept
      : x(initPos[0]), y(initPos[1]), z(initPos[2]) {}
  constexpr position() noexcept {}

  // Operator overloads
  constexpr position operator+(const vec3 &rhs) const noexcept {
    /**
    Addition overload for addition of two positions
    */
    return position(x + rhs[0], y + rhs[1], z + rhs[2]);
  }

  constexpr position operator-(const vec3 &rhs) const noexcept {
    /**
    Subtraction overload for subtraction of two positions
    */
    return position(x - rhs[0], y - rhs[1], z - rhs[2]);
  }

  // Implicit type conversions
  constexpr operator std::array<double, 3>() const noexcept {
    return {x, y, z};
  }

  // Member functions
  constexpr position multiply(double rhs) const noexcept {
    /**
    Multiplication by a constant
    */
    return position(x * rhs, y * rhs, z * rhs);
  }

  constexpr position divide(double rhs) const noexcept {
    /**
    Division by a constant
    */
    return position(x / rhs, y / rhs, z / rhs);
  }

  constexpr position dot(const vec3 &argPos) const noexcept {
    /**
    Performs the vector dot product with the argument
    */
    return position(x * argPos[0], y * argPos[1], z * argPos[2]);
  }
};

/** \struct momentofinertia
\brief Ixx, Iyy, Izz moments
second moment of inertia for the object in three axes
*/
struct momentofinertia {
  // Class variable initialisers
  double Ixx;
  double Iyy;
  double Izz;

  // Member functions
  /** \fn addMass(I, m, r)
  Adds moment of inertia of a component with a mass m at position r from the cg
  */
  constexpr void addMass(const momentofinertia &I, double m,
                         const position &r) noexcept {

    Ixx += I.Ixx + m * r.y * r.y * r.z * r.z;
    Iyy += I.Iyy + m * r.z * r.z * r.x * r.x;
    Izz += I.Izz + m * r.x * r.x * r.x * r.x;
  }
};

/** \struct quaternion
\brief quaternion format
quaternion struct for quaternion rotations
*/
struct quaternion {
  // Class variable initialisers
  double qw = 0;
  double qx = 1;
  double qy = 0;
  double qz = 0;

  // Initialisers
  constexpr quaternion() noexcept {}

  constexpr quaternion(double initW, double initX, double initY,
                       double initZ) noexcept
      : qw(initW), qx(initX), qy(initY), qz(initZ) {}

  /// initialises a quaternion of an equivalent vector and angle

  constexpr quaternion(const vec3 &axis, double angle) noexcept {
    qw = cos(angle * 0.5);
    qx = sin(angle * 0.5) * axis[0];
    qy = sin(angle * 0.5) * axis[1];
    qz = sin(angle * 0.5) * axis[2];
  }
  /// initialises a quaternion of the rotation between two vectors
  constexpr quaternion(const vec3 &arg1, const vec3 &arg2) noexcept {
    vec3 cross = arg1.cross(arg2);
    qw = sqrt(arg1.mag() * arg2.mag()) + arg1.dot(arg2);
    qx = cross[0];
    qy = cross[1];
    qz = cross[2];
  }

  // Member functions
  /** \fn conjugate()
  Finds the conjugate of the quaternion
  */
  constexpr quaternion conjugate() const noexcept {
    return quaternion(qw, -qx, -qy, -qz);
  }

  /** \fn hprod(arg)
  Performs the hamiltonion product of this and the argument quaternion 'arg'
  */
  constexpr quaternion hprod(const quaternion &arg) const noexcept {

    vec3 thisU = vec3(qx, qy, qz); ///< The ijk vector of this quaternion
    vec3 argU =
        vec3(arg.qx, arg.qy, arg.qz); ///< The ijk vector of the quaternion arg

    vec3 qwDotArg = argU * qw;
    vec3 argwDotThis = thisU * arg.qw;
    vec3 thisArgCross = thisU.cross(argU);

    return quaternion(
        qw * arg.qw - thisU.dot(argU),
        qwDotArg[0] + argwDotThis[0] + thisArgCross[0],
        qwDotArg[1] + argwDotThis[1] + thisArgCross[1],
        qwDotArg[2] + argwDotThis[2] +
            thisArgCross
                [2]); // https://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf
  }

  /** \fn toMat()
  converts to Matrix
  */
  constexpr std::array<std::array<double, 3>, 3> toMat() const noexcept {
    /*
    Information for conversion available at:
    https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    */
    std::array<std::array<double, 3>, 3> result = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    result[0] = {{1 - 2 * qy * qy - 2 * qz * qz, 2 * qx * qy - 2 * qz * qw,
                  2 * qx * qz + 2 * qy * qw}};
    result[1] = {{2 * qx * qy + 2 * qz * qw, 1 - 2 * qx * qx - 2 * qz * qz,
                  2 * qy * qz - 2 * qx * qw}};
    result[2] = {{2 * qx * qz - 2 * qy * qw, 2 * qy * qz + 2 * qx * qw,
                  1 - 2 * qx * qx - 2 * qy * qy}};
    return result;
  }
  /** \fn quaternionDerivative(argQuat, argRate)
  Takes the derivative of an argument quaternion and corresponding argument rate
  */
  constexpr quaternion
  quaternionDerivative(quaternion argQuat, const vec3 &argRate) const noexcept {
    quaternion dotQuat;

    argQuat.qw =
        sqrt(1 - pow2(argQuat.qx) - pow2(argQuat.qy) - pow2(argQuat.qz));

    dotQuat.qx = 0.5 * (argQuat.qw * argRate[0] - argQuat.qz * argRate[1] +
                        argQuat.qy * argRate[2]);

    dotQuat.qy = 0.5 * (argQuat.qz * argRate[0] + argQuat.qw * argRate[1] -
                        argQuat.qx * argRate[2]);

    dotQuat.qz = 0.5 * (-argQuat.qy * argRate[0] + argQuat.qx * argRate[1] -
                        argQuat.qw * argRate[2]);

    return dotQuat;
  }

  /** \fn rotate(arg)
  Rotates an argument vector3 by this quaternion
  */
  constexpr vec3 rotate(const vec3 &arg) const noexcept {
    /*
    Rotates vec3 arg using by this quaternion
    */
    quaternion argq = quaternion(0, arg[0], arg[1], arg[2]);

    quaternion returnq = hprod(argq.hprod(
        conjugate())); // for this quaternion p and argument quaternion q, this
                       // is the quaternion multiplication pqp*

    return vec3(returnq.qx, returnq.qy, returnq.qz);
  }

  constexpr std::array<double, 4> getAxisAngle() const noexcept {

    vec3 axis = vec3(qx, qy, qz) / vec3(qx, qy, qx).mag();
    return {acos(qw) * 2, axis[0], axis[1], axis[2]};
  }
};

/** \struct ftModel
\brief force-torque model
Represents the vector of forces and moments of an actuator at maximum actuation
*/
struct ftModel {

  /// Generic initialiser
  std::array<double, 6> ftVec{0, 0, 0, 0, 0, 0};

  /// Initialisers
  constexpr ftModel(double Fx, double Fy, double Fz, double Txx, double Tyy,
                    double Tzz) noexcept
      : ftVec({Fx, Fy, Fz, Txx, Tyy, Tzz}) {
    /*
    Generic initialiser for the vector, directly assigning each component
    */
  }

  /** initialiser for a thruster with thrust magnitude and direction of action
   */
  constexpr ftModel(const vec3 &maxThrust, const vec3 &thrusterPos) noexcept {

    ftVec[0] =
        maxThrust[0]; // Resultant force on the object cg is off centre force
    ftVec[1] = maxThrust[1];
    ftVec[2] = maxThrust[2];

    vec3 torque = thrusterPos.cross(maxThrust);

    ftVec[3] =
        torque[0]; ///< Resultant force on the object cg is off centre force
    ftVec[4] =
        torque[1]; ///< Resultant force on the object cg is off centre force
    ftVec[5] =
        torque[2]; ///< Resultant force on the object cg is off centre force
  }

  constexpr ftModel(const std::array<double, 6> &initFTV) noexcept
      : ftVec(initFTV) {}

  // Operators
  constexpr double operator[](size_t i) const noexcept { return ftVec[i]; }
  /**
  Addition overload for addition of two ftModels
  */
  constexpr ftModel operator+(const ftModel &rhs) const noexcept {
    return ftModel(ftVec[0] + rhs[0], ftVec[1] + rhs[1], ftVec[2] + rhs[2],
                   ftVec[3] + rhs[3], ftVec[4] + rhs[4], ftVec[5] + rhs[5]);
  }
  /**
  Subtraction overload for addition of two ftModels
  */
  constexpr ftModel operator-(const ftModel &rhs) const noexcept {

    return ftModel(ftVec[0] - rhs[0], ftVec[1] - rhs[1], ftVec[2] - rhs[2],
                   ftVec[3] - rhs[3], ftVec[4] - rhs[4], ftVec[5] - rhs[5]);
  }
  /**
  Multiplication overload for addition of two ftModels
  */
  constexpr ftModel operator*(double rhs) const noexcept {

    return ftModel(ftVec[0] * rhs, ftVec[1] * rhs, ftVec[2] * rhs,
                   ftVec[3] * rhs, ftVec[4] * rhs, ftVec[5] * rhs);
  }

  // Implicit type converters
  /**
  Adds support for implicit conversion from ftModel to std::array<double, 6>
  */
  constexpr operator std::array<double, 6>() const noexcept { return ftVec; }

  // Member functions
  constexpr double mag() const noexcept {
    return sqrt(pow2(ftVec[0]) + pow2(ftVec[1]) + pow2(ftVec[2]) +
                pow2(ftVec[3]) + pow2(ftVec[4]) + pow2(ftVec[5]));
  }
  /** \fn normalise()
  normalises the ftModel*/
  constexpr ftModel normalise() const noexcept {
    double invmag = 1 / mag();
    return std::array<double, 6>{ftVec[0] * invmag, ftVec[1] * invmag,
                                 ftVec[2] * invmag, ftVec[3] * invmag,
                                 ftVec[4] * invmag, ftVec[5] * invmag};
  }
  /** \fn getDominantAxis()
  returns the dominant axis of the thrust*/
  constexpr int getDominantAxis() const noexcept {
    double maxVal = 0;
    int maxValAxis = 0;
    for (int i = 0; i < 6; i++) {
      if (abs(ftVec[i]) > maxVal) {
        maxVal = abs(ftVec[i]);
        maxValAxis = i;
      }
    }
    return maxValAxis;
  }
};

/** \struct rotStates
\brief rotation states
a struct composed of 2 vec3s */
struct rotStates {
  /// @param omega body rates
  vec3 omega; // body rates
  /// @param q the vector part of the quaternion
  vec3 q; // quaternion vector part
};

/** \struct posStates
\brief position states
a struct composed of the earth-centered inertial data */
struct posStates {
  /// @param r earth-centered inertial position
  vec3 r; // eci position (m)
  /// @param v earth-centered inertial velocity
  vec3 v; // eci velocity (m/s)
  /// @param m mass of spacecraft
  double m; // spacecraft mass (kg)
};

/** \struct orbParam
\brief orbital parameters
structure of the orbital parameters */
struct orbParam {
  /// @param sma semi major axis (m)
  long int sma;
  /// @param ecc eccentricity
  double ecc;
  /// @param inc inclination (rad)
  double inc;
  /// @param asc argument of periapsis (rad)
  double asc;
  /// @param aop longitude of the ascending node (rad)
  double aop;
  // double meanAnom; // mean anomaly (rad)
  // double eccAnom; // eccentric anomaly (rad)
  /// @param truAnom true anomaly (rad)
  double truAnom;

  // member functions
  /** \fn meanToTrue(meanAnomaly)
  @param[in] meanAnomaly input mean anomaly
  Converts mean anomaly to true anomaly
  */
  double meanToTrue(double meanAnomaly){
    // converts mean anomaly to true anomaly, used for calculating positions at
    // different times
    double dTheta = 0;
    if (ecc < 1.0) { // newton raphson's method must be used for higher
                          // eccentricities, e>1 is a parabolic orbit

      double Beta = ecc/(1+sqrt(1+pow2(ecc)));
      double EccAnomaly = meanAnomaly + ((ecc * sin(meanAnomaly)) / ((cos(ecc) - ((M_PI_2 - ecc) * sin(ecc))) + meanAnomaly * sin(ecc)));
      double dE = 10e10;

      while (abs(dE) > 10e-10) {
        dE = (EccAnomaly - ecc * sin(EccAnomaly) - meanAnomaly) / (1 - ecc * cos(EccAnomaly));
        EccAnomaly -= dE;
      };

      dTheta = atan2(sqrt(1-pow2(ecc))*sin(EccAnomaly), cos(EccAnomaly)-ecc);
    };
    
    return dTheta;
  };

  /** \fn trueToMean()
  converts the true anomaly to mean anomaly
  */
  constexpr double trueToMean(double deltaT) noexcept {
    double meanAnomaly = 0;
    double a = sqrt(1 - pow2(ecc)) * sin(deltaT);
    double b = 1 + ecc * cos(deltaT);

    meanAnomaly = atan2((a / b), ((ecc + cos(deltaT)) / b)) - ecc * (a / b);

    return meanAnomaly;
  };

  double MeanOrbitalMotion(){
    return sqrt(mu / pow3(sma)); // unfuck this later when you learn how to code
                                 // not like an ape
  }

  double ArcTimeCalc(double Theta){
     return (1 / pow2(1 + ecc * cos(Theta)));
  }

  // constexpr double RK4(double dTheta, double stepSize) noexcept {
  //   double k1, k2, k3, k4, y0, y1, y2, y3, y4, yn, tn, y = 0;
  //   y0 = ArcTime(truAnom);
  //   for (double x = truAnom; x <= truAnom + dTheta; x += stepSize) {
  //     k1 = ArcTimeDerivative(x);
  //     y1 = y0 + k1 * stepSize / 2;
  //     k2 = ArcTimeDerivative(x + stepSize / 2);
  //     y2 = y0 + k2 * stepSize / 2;
  //     k3 = ArcTime(x + stepSize / 2);
  //     y3 = y0 + k3 * stepSize;
  //     k4 = ArcTime(x + stepSize);
  //     yn = y0 + stepSize * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  //     y0 = yn;
  //   }
  //   return yn;
  // }

  double CompositeTrapezoid(double dTheta){//fucked somehow idk
    double ArcTime;
    double n = 10000; //subinterval count
    double h = dTheta/n;
    double f_a = ArcTimeCalc(truAnom)/2;
    double f_b = ArcTimeCalc(truAnom+dTheta)/2;
    double f_k=0;
    for (int k=1; k<n; k++){
      f_k+=ArcTimeCalc(truAnom+k*h);
    }

    ArcTime = h * (f_a + f_k + f_b)/(pow2(mu)/pow3(AngularMomentum()));
    return ArcTime;
  }

  double AngularMomentum(){
    double h = sqrt(mu*sma*(1-pow2(ecc)));
    return h;
  }

  constexpr double KeplersEqn(double dMean) const noexcept {
    return dMean / sqrt(mu / pow3(sma));
  }
};

/** \struct pcs
\brief PCS coordinate system
Orbital position of satellite in Perifocal Co-Ordinate System - satellite
centred
*/
struct pcs {
  // p - points towards periapsis of orbit
  // q - right hand angle along orbital plane
  // w - normal to orbital plane
  /// @param rPCS position in perifocal coordinate system
  vec3 rPCS;
  /// @param vPCS velocity in perifocal coordinate system
  vec3 vPCS;
};

/** \struct eci
\brief ECI coordinate system
Orbital position of satellite in Earth Centred Inertial Co-Ordinate System
*/
struct eci {
  // i - vector from Earth to sun on J2000, 2000-01-01 at 12:00 TT
  // j - orthogonal towards Equator
  // k - passes through Celestial North Pole
  /// @param rIJK position in ECI coordinates
  vec3 rIJK;
  /// @param vPCS velocity in ECI coordinates
  vec3 vIJK;
};

/** \struct ecef
\brief ECEF coordinate system
Orbital position of satellite in Earth Centred Earth Fixed Co-Ordinate System
*/
struct ecef {
  // x - vector passing through Greenwich Meridian
  // y - orthogonal towards Equator
  // z - passes through Celestial North Pole
  /// @param rXYZ position in ECEF coordinates
  vec3 rXYZ;
  /// @param vXYZ velocity in ECEF coordinates
  vec3 vXYZ;
};

/** \struct ned
\brief NED coordinate system
Orbital position of satellite in North East Down Co-Ordinate System - satellite
centred
*/
struct ned {
  // n - points North
  // e - points East
  // d - points down
  /// @param rNED position in NED coordinates
  vec3 rNED;
  /// @param vNED velocity in NED coordinates
  vec3 vNED;
};

/** \struct enu
\brief ENU coordinate system
Orbital position of satellite in East North Up Co-Ordinate System - Earth
centred
*/
struct enu {
  // e - points East
  // n - points North
  // u - points Up
  /// @param rENU position in ENU coodinates
  vec3 rENU;
  /// @param vENU velocity in ENU coordinates
  vec3 vENU;
};

/** \struct thcs
\brief THCS coordinate system
Orbital position of satellite in Topocentric Horizon Co-Ordinate System -
centred on ground point
*/
struct thcs {
  // s - points South
  // e - points East
  // z - points up
  /// @param rSEZ position in SEZ coordinates
  vec3 rSEZ;
  /// @param vSEZ velocity in SEZ coordinates
  vec3 vSEZ;
};

/** \struct lla
\brief LLA coordinate system
Ground Sub-Vehicle Point on Earth's surface
*/
struct lla {
  /// @param lat geocentric longitude
  double lat;
  /// @param lon geocentric latitude
  double lon;
  /// @param alt height above WGS84 ellipsoid
  double alt;
};

/** \struct ear
\brief EAR coordinate system
Orbital position of satellite in Elevation Azimuth Range Co-Ordinate System -
centred on ground point
*/
struct ear {
  /// @param e elevation (rad)
  double e;
  /// @param a azimuth (rad)
  double a;
  /// @param r points upwards
  double r;
};

/** \struct vnb
\brief VNB coordinate system
Orbital velocity of satellite in Velocity Normal Bi-Normal axes - used for
orbital maneuvers
*/
struct vnb {
  // v - velocity vector - prograde
  // n - normal vector - normal to orbital plane
  // b - bi-normal vector - orthangonal to x and y, in the direction of the
  // angular momentum vector
  /// @param vVNB velocity in VNB coordinate system
  vec3 vVNB;
};

/** \struct powermodel
\brief power states structure
Structure to hold different power states of a component
*/
struct powermodel {
  // 4 different states power states
  /// @param off power consumption during off state
  double off = 0;
  /// @param idle power consumption during idle state
  double idle = 0;
  /// @param use power consumption during use state
  double use = 0;
  /// @param max maxmium power consuption of component
  double max = 0;
};
} // namespace osc

#endif // OSCTYPES_H
#ifndef PLANET_H
#define PLANET_H

#include <vector>
#include <array>
#include <math.h>

//default values for Earth
struct celestial{
    /// @param sgp standard gravitational parameter mu
    double sgp;
    /// @param sMa equatorial radius of planet
    double sMa; 
    /// @param flat ellipticity of planet's surface
    double flat; 
    /// @param sma polar radius of planet
    double sma; 
    /// @param ecc ecccentricity of planet's surface
    double ecc; 
    /// @param J2 second zonal gravitational harmonic
    double J2; 
    celestial(double initsgp, double initsMa, double initfla, double initsma, double initecc, double initJ2):
        sgp(initsgp), sMa(initsMa), flat(initfla), sma(initsma), ecc(initecc), J2(initJ2){};
};

// celestial earth(3.986004418e14, 6378137, 1/298.257223563,   6356752.314,    0.08181919084,  1082.63e-6);
// celestial venus(3.24859e14,     6051800, 0,                 6051800,        0,              4.458e-6);
// celestial mars(4.282837e13,     3396200, 0.00589,           3376200,        0.1083757717,   1960.45e-6);
// celestial moon(4.9048695e12,    1738100, 0.0012,            1736000,        0.04897509571,  202.7e-6);

const celestial planet(3.986004418e14, 6378137, 1/298.257222101,   6356752.314,    0.08181919084,  1082.63e-6); // defined as earth

#endif // PLANET_H
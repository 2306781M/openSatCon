#include <iostream>

#include "osctypes.hpp"

namespace osc{
double g                =   9.80665;
int h_f                 =   283000;         //final height of ASAT-M
int t_f                 =   168;            //time to reach final height
int V_f                 =   3400;           //final velocity of ASAT-M
int Mo                  =   19000;          //initial mass of ASAT-M
double Lambda           =   19/(19-16.7);   //Payload fraction values 
int Isp                 =   280;            //standard value for solid fuel motor
double ModelAccuracy    =   0.001;          //modelling value
double mu               =   3.986012e14;    //keeping large values to avoid
int EarthRadius         =   6371000;        //unit changing mid-model
int AltAtmosphere       =   84852+EarthRadius;

//needs to create a KOE object with these parameters
    // Sat_Alt=250000; Sat_Ecc=0;      Sat_Inc=0;
    // Sat_RAAN=0;     Sat_ArgPer=0;   Sat_MeaAno=0;
    // Sat_SMA=EarthRadius+Sat_Alt;
    // SatKOE=[Sat_SMA Sat_Ecc Sat_Inc Sat_RAAN Sat_ArgPer Sat_MeaAno];

//defender satellite initial keplerian orbital elements
    // CaSat_Alt_i=250000; CaSat_Ecc_i=0;      CaSat_Inc_i=0;
    // CaSat_RAAN_i=0;     CaSat_ArgPer_i=0;   CaSat_MeaAno_i=0;
    // CaSat_SMA_i=EarthRadius+CaSat_Alt_i;
    // CaSatKOEi=[CaSat_SMA_i CaSat_Ecc_i CaSat_Inc_i CaSat_RAAN_i CaSat_ArgPer_i CaSat_MeaAno_i];

//ASAT-M Model Values
double mo=((Mo/exp(((V_f/g)+t_f)/280))-Mo)/(-t_f);     //fuel mass flow Rate
double p_f=1-((Mo-mo*t_f)/Mo)*(log(Mo/(Mo-mo*t_f))+1); //range function p(f)
double TWR=(g*pow2(Isp)*p_f)/(h_f+0.5*g*pow2(t_f));            //thrust to weight Ratio

//ASAT-M Modelling with time
//t=ModelAccuracy:ModelAccuracy:t_f;              //sets up matrix for all times
int t;
 
double p=1-((Mo-mo*t)/Mo)*(log(Mo/(Mo-mo*t))+1);  //range function p(f) again
double h=((g*pow2(Isp)/TWR)*p-(0.5*g*pow2(t)));              //altitude as function of time
double V=g*(Isp*log(Mo/(Mo-mo*t))+1);                 //velocity as function of time
double h_c=h+(pow2(V))/(2*g);      //culmination alt as func of time

//for t_burnout=t_f+ModelAccuracy:ModelAccuracy:1000 //gives height values after motor burnout
double t_burnout;
//double h = h(t_f/ModelAccuracy)+V(t_f/ModelAccuracy)*(t_burnout-t_f)-(g/2*(t_burnout-t_f));//altitude after burnout


//HRCT = find(h>Sat_Alt, 1);                      //first time index where ASAT-M is at target satellite
//LRCT = (find(h_c>100000, 1)-1);                 //first time index where culmination height is above 100km
//ASAT_AltAtLRCT=h(LRCT);                         //Altitude of ASAT-M at Low Risk Critical Time
}
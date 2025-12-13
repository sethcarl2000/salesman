
#include "simulated_annealing.hpp"
#include "CityCoord.hpp"
#include <cmath> 

namespace {
  constexpr double deg_to_rad = 0.01745329251; //pi/180
}

//__________________________________________________________________________________
// given two cities, return the great-circle min-distance between them, in radians. 
double CityDistance(const CityCoord& c1, const CityCoord& c2)
{
  using namespace std; 

  //difference betweem the two phi's
  const double cos_phi1_2 = cos(deg_to_rad*(c1.lon - c2.lon));
  
  //thetas (with 0 being the north pole)
  const double sin_th1 = sin(deg_to_rad*(90. - c1.lat)); 
  const double sin_th2 = sin(deg_to_rad*(90. - c2.lat)); 

  return acos(
    (cos_phi1_2 * sin_th1 * sin_th2) + sqrt((1. - sin_th1*sin_th1)*(1. - sin_th2*sin_th2))
  );
} 

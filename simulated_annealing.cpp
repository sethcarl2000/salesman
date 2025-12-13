
#include "simulated_annealing.hpp"
#include "CityCoord.hpp"
#include <cmath> 
#include <stdexcept> 
#include <sstream> 
#include <random> 
#include <functional> 

namespace {
  constexpr double pi = 3.14159265359; 
  constexpr double deg_to_rad = pi/180.; 
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
//_________________________________________________________________________________


//_________________________________________________________________________________
// Create an array of 'fake' cities, with a given generation pattern. 
// 
//  Valid generation keywords are: 
//    -   "random"      -cities are randomly distributed over the earth's surface
//    -   "equatorial"  -cities are randomly generated along the equator
//    -   "polar"       -cities are randomly generated along a great circle that intersects with the north pole
std::vector<CityCoord> MakeFakeCities(const size_t N, const std::string type="random", std::random_device* rd=nullptr)
{
  using namespace std; 

  //check if the argument provided is valid or not
  bool valid_arg{false}; 

  for (const string& arg : {"random", "equatorial", "polar"}) { 
    if (arg==type) { valid_arg=true; break; } 
  }

  //throw an exception and quit if the generation type is invalid
  if (!valid_arg) {   
    ostringstream oss; oss << "in <" << __func__ << ">: Invalid argument passed for generation type: '" << type << "'. See fcn body for valid arguments."; 
    throw invalid_argument(oss.str()); 
    return {}; 
  }

  //if we weren't passed a random_device, create our own.
  // (I'm gonna leave this ptr dangling, which is bad practice, but I can't imagine that this fcn will
  //  be invoked more than once per execution... so it's probably fine.) 
  if (!rd) rd = new random_device; 

  vector<CityCoord> cities; cities.reserve(N); 

  //we will use this to pick a random latitude/longitude
  if (type=="random") {
    
    //this will give us random, gauss-distributed numbers with mean=0 and sigma=1 
    auto gaus = std::bind( normal_distribution{0., 1.}, mt19937((*rd)()) ); 

    for (size_t c=0; c<N; c++) {
  
      //Create a randomly-oriented unit-vector 
      double norm =0.; 
      double X[3];
      for (int i=0; i<3; i++) {
        X[i] = gaus(); 
        norm += X[i]*X[i]; 
      }
      norm = 1./sqrt(norm); 
      for (int i=0; i<3; i++) X[i] *= norm; 

      //Now, we can convert this unit vector into a lat/lon
      double cos_theta = X[2]; 
      double cos_phi   = X[0] / sqrt(1. - cos_theta*cos_theta); 

      double lat = (pi/2. - acos(cos_theta)) / deg_to_rad; 
      double lon = acos(cos_phi) * (X[1]/fabs(X[1])) / deg_to_rad;
      
      //add it to the list of randomly-generated cities
      cities.push_back(CityCoord{.lat=lat, .lon=lon}); 
    }
    return cities; 
  } 
  
  

  return {}; 
}



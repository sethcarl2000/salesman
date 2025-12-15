#ifndef simulated_annealing_HPP
#define simulated_annealing_HPP

#include "CityCoord.hpp"
#include <string> 
#include <vector> 
#include <random> 

//__________________________________________________________________________________
// given two cities, return the great-circle min-distance between them, in radians. 
double CityDistance(const CityCoord& c1, const CityCoord& c2); 

//__________________________________________________________________________________
// Create an array of 'fake' cities, with a given generation pattern. 
std::vector<CityCoord> MakeFakeCities(const size_t N, const std::string type="random", std::random_device* rd=nullptr); 




#endif 
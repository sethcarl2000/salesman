
#include "simulated_annealing.hpp"
#include "CityCoord.hpp"
#include <cmath> 
#include <stdexcept> 
#include <sstream> 
#include <random> 
#include <functional> 
#include <algorithm> 
#include <iostream> 

namespace {
  constexpr double pi = 3.14159265359; 
  constexpr double deg_to_rad = pi/180.; 
  constexpr double half_pi = pi/2.; 
}

//__________________________________________________________________________________
// given two cities, return the great-circle min-distance between them, in radians. 
double CityDistance(const CityCoord& c1, const CityCoord& c2)
{
  using namespace std; 

  //difference betweem the two phi's
  const double cos_phi1_2 = cos(c1.lon - c2.lon);
  
  //thetas (with 0 being the north pole; same as for spherical coordinates)
  const double sin_th1 = sin(half_pi - c1.lat); 
  const double sin_th2 = sin(half_pi - c2.lat); 

  const double cos_th1 = cos(half_pi - c1.lat); 
  const double cos_th2 = cos(half_pi - c2.lat); 
  
  return acos((cos_phi1_2 * sin_th1 * sin_th2) +  (cos_th1 * cos_th2)); 

} 
//_________________________________________________________________________________


//_________________________________________________________________________________
// Create an array of 'fake' cities, with a given generation pattern. 
// 
//  Valid generation keywords are: 
//    -   "random"      -cities are randomly distributed over the earth's surface
//    -   "equatorial"  -cities are randomly generated along the equator
std::vector<CityCoord> MakeFakeCities(const size_t N, const std::string type)
{
  using namespace std; 

  //check if the argument provided is valid or not
  bool valid_arg{false}; 

  for (const string& arg : {"random", "equatorial"}) { 
    if (arg==type) { valid_arg=true; break; } 
  }

  //throw an exception and quit if the generation type is invalid
  if (!valid_arg) {   
    ostringstream oss; oss << "in <" << __func__ << ">: Invalid argument passed for generation type: '" << type << "'. See fcn body for valid arguments."; 
    throw invalid_argument(oss.str()); 
    return {}; 
  }

  random_device rd; 

  vector<CityCoord> cities; cities.reserve(N); 

  //we will use this to pick a random latitude/longitude
  if (type=="random") {
    
    //this will give us random, gauss-distributed numbers with mean=0 and sigma=1 
    auto gaus = std::bind( normal_distribution{0., 1.}, mt19937(rd()) ); 

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

  //randomly distributed cities along the equator
  if (type=="equatorial") {

    //randomly generates longitude
    auto rand_longitude = std::bind( uniform_real_distribution<double>{-180., +180.}, mt19937(rd()) ); 

    for (size_t i=0; i<N; i++) {
      cities.push_back(CityCoord{.lat=0., .lon=rand_longitude()}); 
    }
    return cities; 
  }

  return {}; 
}

//__________________________________________________________________________________
// Sort this list of cities, starting with the first city, by finding the next-closest city for each neighbor. 
std::vector<int> SortCities(const std::vector<CityCoord>& cities)
{
  using namespace std; 

  const int n_cities = cities.size(); 

  //create the default ordering of cities
  vector<int> rand_order, order; 
  rand_order.reserve(n_cities); 
  order     .reserve(n_cities); 

  for (int i=0; i<n_cities; i++) rand_order.push_back(i); 

  std::random_device rd; 

  std::shuffle( order.begin(), order.end(), mt19937(rd()) );

  //now, let's sort the cities by their nearest-neighbors

  //take one element from the 'rand_order' vector, and erase it from the other. 
  auto insert_element = [&rand_order,&order](int index){
    
    if (index >= (int)rand_order.size() || index < 0) {
      throw logic_error("invalid index passed.");
      return; 
    }
    order.push_back(rand_order[index]); 
    rand_order.erase(rand_order.begin()+index); 
    return; 
  };

  //take the first element of 'rand_order' as our starting city
  insert_element(0); 

  while (rand_order.size() > 0) {

    const auto& city = cities[order.back()]; 

    double min_dist = 1.e30;
    int min_city_index = -1; 

    for (size_t i=0; i<rand_order.size(); i++) {

      CityCoord new_city = cities[rand_order[i]];
      
      double dist = CityDistance(new_city, city); 
    
      if (dist < min_dist) {
        min_city_index = i; 
        min_dist = dist; 
      }
    }

    insert_element(min_city_index); 
  }

  return order; 
}

//__________________________________________________________________________________
Annealer::Annealer(const std::vector<CityCoord>& _c) 
  : fCities{_c}, n_cities{(int)_c.size()} 
{
  using namespace std;
  
  //initialize the random-number generators
  random_device rd; 
  fGenerator = mt19937(rd());   
  fUniform_dist = uniform_real_distribution<double>(0., 1.); 
  fRand_index_dist = uniform_int_distribution<int>(0., n_cities-1); 

  fOrder = SortCities(fCities); 

  //compute distances of all cities 
  //yes, this array is enormous. 
  fDistances = new double[n_cities*n_cities]; 
  for (int i1=0; i1<n_cities; i1++) {
    for (int i2=0; i2<n_cities; i2++) {

      fDistances[ i1*n_cities + i2 ] = CityDistance( fCities[i1], fCities[i2] ); 
    }
  }
}; 

//___________________________________________________________________________________
//randomly shuffle order of cities
void Annealer::Shuffle() 
{
  std::shuffle( fOrder.begin(), fOrder.end(), fGenerator ); 
}

//___________________________________________________________________________________
//get total length of current ordering (in radians)
double Annealer::Get_total_length() const
{
  double length = 0.; 
  for (size_t i=0; i<fCities.size()-1; i++) {

    length += CityDistance_arr( fOrder[i], fOrder[i+1] );
  }
  length += CityDistance_arr( fOrder.back(), fOrder.front() );

  return length;
}

//__________________________________________________________________________________
// annealing update where the proposed update is a swap of the index ordering 
double Annealer::Swap_update(double T, const unsigned long n_swaps) 
{
  using namespace std; 
  //check the 'energy change' (change in total angle) between the current state and the proposed update. 

  //number of proposed swaps that are accepted. 
  unsigned long n_accepted=0;

  //loop through all the swaps 
  for (unsigned long i=0; i<n_swaps; i++) {

    //get two (distinct!) random indices
    int ind1, ind2; 

    ind1 = rand_index(); 
    do { ind2 = rand_index(); } while (ind1==ind2);
    
    //check to see if the indices are adjacent
    int sep = abs( ind1 - ind2 ); 
    bool is_adjacent = (sep==1 || sep==n_cities-1); 

    //this will be the change in length (in radians). 
    // d_energy < 0 means the suggested update is a decrease in length,
    // d_energy > 0 means the suggested update is an increase in length.  
    double d_energy; 

    if (is_adjacent) { //these two cities are currently adjacent in our ordering of cities

      int i_lo = min<int>(ind1, ind2); 
      int i_hi = max<int>(ind1, ind2); 

      //first compute the original length.
      //we will cheat here; if we swap the order of two adjacent cities, then the __distance between them__ does not change; 
      // the only change in length comes from swapping their order: 
      /*  
            ---- C0   C2   
                  |  / |
                  | /  |
                 C1   C3 ----- 

      swapping the order of C1 & C2: 

            ---- C0---C2   
                    / 
                  /  
                C1---C3 ----- 
      */
      d_energy = 
        CityDistance_arr( 
          i_lo-1>=0       ? fOrder[i_lo-1] : fOrder.back(), 
          fOrder[i_hi] 
        ) + 
        CityDistance_arr( 
          i_hi+1<n_cities ? fOrder[i_hi+1] : fOrder.front(), 
          fOrder[i_lo] 
        ); 

      //now, compute what these distances are in the original (current) ordering. 
      d_energy -= 
        CityDistance_arr( 
          i_lo-1>=0       ? fOrder[i_lo-1] : fOrder.back(), 
          fOrder[i_lo] 
        ) + 
        CityDistance_arr( 
          i_hi+1<n_cities ? fOrder[i_hi+1] : fOrder.front(), 
          fOrder[i_hi] 
        ); 
      

    } else { //this is for cities that are NOT currently adjacent. 
        
      //this only works for NON-ADJACENT cities (adjaceny defined by the current ordering of cities, NOT by physical proximity)
      auto check_swapped_distances = [&](int index_old, int index_new)
      {
        return
          CityDistance_arr( 
            index_old-1>=0       ? fOrder[index_old-1] : fOrder[n_cities-1], 
            fOrder[index_new] 
          ) +  
          CityDistance_arr( 
            fOrder[index_new], 
            index_old+1<n_cities ? fOrder[index_old+1] : fOrder[0] 
          ); 
      };

      d_energy = 
        check_swapped_distances(ind1, ind2) + 
        check_swapped_distances(ind2, ind1); 

      d_energy -= 
        check_swapped_distances(ind1, ind1) + 
        check_swapped_distances(ind2, ind2); 
    }

    //now randomly decide if the swap will be accepted
    if (d_energy < 0. || exp( -d_energy/T ) > rand_uniform()) {

      //then, accept the new update. swap the orders of the cities
      int index_new_1 = fOrder[ind2]; 
      int index_new_2 = fOrder[ind1]; 

      fOrder[ind1] = index_new_1;  
      fOrder[ind2] = index_new_2; 

      n_accepted++; //new update accepted
    } 
  }

  return ((double)n_accepted)/((double)n_swaps); 
} 



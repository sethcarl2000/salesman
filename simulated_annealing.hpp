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
std::vector<CityCoord> MakeFakeCities(const size_t N, const std::string type="random"); 

//__________________________________________________________________________________
// Sort this list of cities, starting with the first city, by finding the next-closest city for each neighbor. 
std::vector<int> SortCities(const std::vector<CityCoord>& cities);

//is this a verb? Anneal-ing-er... 
class Annealer {
public: 
    Annealer(const std::vector<CityCoord>& _c);

    //update that considers a swap between two randomly-chosen cities
    bool Swap_update(double beta); 

    //update that picks two random cities, and considers putting them next-to-each other as the update. more expensive! 
    bool Pinch_update(double beta); 


    //Get the total length of the current path, in radians
    double Get_total_length() const; 

    //get the ordering of cities
    std::vector<int> Get_city_order() const { return fOrder; }

    //get the list of cities
    std::vector<CityCoord> Get_cities() const { return fCities; }
    
private: 

    //compute the length of a walk between several cities, ordered by the indices provided.
    double Compute_chain_length(const std::vector<int>& indices) const; 

    std::vector<CityCoord> fCities;
    const int n_cities; 
    std::vector<int> fOrder; 
    std::uniform_int_distribution<int> fRand_index_dist; 
    std::uniform_real_distribution<double> fUniform_dist; 
    std::mt19937 fGenerator; 

    //random number generators 
    inline int rand_index() { return fRand_index_dist(fGenerator); }
    inline double rand_uniform() { return fUniform_dist(fGenerator); }
};


#endif 
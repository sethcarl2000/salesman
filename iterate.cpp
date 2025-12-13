// example code to read in a data file of city lat,long coordinates

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random> 
#include "simulated_annealing.hpp"
#include "microsecond_execution_timer.hpp"
#include "CityCoord.hpp"
#include <cstdio> 
#include <iostream> 

using namespace std; 

// simple structure to store city coordinates
// could also use std::pair<double,double> 
// or define a class

//__________________________________________________________________________________________
// fill the array of city locations
vector<CityCoord> GetData(char* fname){
  
  FILE* fp=fopen(fname,"r");
  
  const int bufsiz=1000;
  char line[bufsiz+1];
  int ncity=0;
  
  vector<CityCoord> cities; 

  while(1){
    CityCoord city; 
    
    fgets(line,bufsiz,fp);
    if (line[0]=='#') continue;  // skip comments
    if (feof(fp)) break;
    // we only scan for two numbers at start of each line
    sscanf(line, "%lf %lf", &city.lon, &city.lat);

    cities.push_back(city); 
  }
  fclose(fp);
  
  return cities;
}
//__________________________________________________________________________________________

int main(int argc, char *argv[]){
  
  const int NMAX=2500;
  
  if (argc<2){
    printf("Please provide a data file path as argument\n");
    return 1;
  }

  auto cities = GetData(argv[1]);

  printf("Read %zi cities from data file\n", cities.size());

  //printf("Longitude  Latitude\n");
  //for (const auto& city : cities) printf("%lf %lf\n",	city.lon, city.lat); 

  long int n_tests = 1e7; 

  random_device rd; 
  mt19937 randgen(rd()); 
  uniform_int_distribution<int> rand_index_generator(0, cities.size()-1); 

  auto rand_index = std::bind( rand_index_generator, std::ref(randgen) ); 

  auto time_with_measure = [&cities, &rand_index, n_tests](){
    for (long int i=0; i<n_tests; i++) {

      const auto& c1 = cities[rand_index()]; 
      const auto& c2 = cities[rand_index()]; 

      double dist = CityDistance(c1, c2); 
    }
  }; 

  auto time_without_measure = [&cities, &rand_index, n_tests](){
    for (long int i=0; i<n_tests; i++) {

      const auto& c1 = cities[rand_index()]; 
      const auto& c2 = cities[rand_index()]; 
    }
  }; 

  printf("measuring execution time with %li events...", n_tests); cout << flush; 

  double microseconds_with_measure    = microsecond_execution_timer(time_with_measure); 
  double microseconds_without_measure = microsecond_execution_timer(time_without_measure); 

  printf(
    "done.\n"
    "Timing of %li events:\n"
    " with measurement    - %.6f us/event\n"
    " without measurement - %.6f us/event\n" 
    "  ~~ so, esitmated time per measurement is:  %.6f us/measurement ~~ \n", 
    
    n_tests, 
    microseconds_with_measure/((double)n_tests),
    microseconds_without_measure/((double)n_tests), 
    (microseconds_with_measure - microseconds_without_measure)/((double)n_tests)
  ); 

  cout << endl; 

  return 0;
}


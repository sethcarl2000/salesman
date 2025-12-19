// example code to read in a data file of city lat,long coordinates

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono> 
#include <random> 
#include <algorithm> 
#include "simulated_annealing.hpp"
#include "microsecond_execution_timer.hpp"
#include "NiceOptionParser.hpp"
#include "CityCoord.hpp"
#include <cstdio> 
#include <iostream> 
#include <iostream> 
#include <TH1D.h> 
#include <TGraph.h> 
#include <TAxis.h>
#include <TCanvas.h> 
#include <TPad.h> 
#include <TStopwatch.h> 
#include <stdexcept> 
#include <fstream> 

using namespace std; 

namespace {
  constexpr double pi = 3.14159265359; 
  constexpr double deg_to_rad = 180./pi; 

  constexpr double earth_radius_km = 6356.7523; 
}
//get the maximum element from a std::vector<T>
template<typename T> T vector_max(const std::vector<T>& V, T max=(T)0.) {
  for (const T& x : V) if (x > max) max = x; 
  return max; 
}

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

    //convert from degrees to radians for the lat/lon 
    city.lon *= deg_to_rad; 
    city.lat *= deg_to_rad; 

    cities.push_back(city); 
  }
  fclose(fp);
  
  return cities;
}
//__________________________________________________________________________________________
//returns the total length of a fully-connected path, in radians
void save_city_csv(const string& path_out, const vector<CityCoord>& cities, const vector<int>& order)
{
  cout << "\nMaking output file '" << path_out << "'..." << flush; 

  //create the output file 
  ofstream outfile(path_out.data(), ios::trunc); 

  outfile << "# longitude latitude\n"; 
  for (int index : order) {
    const auto& city = cities[index];
    
    char line[500]; sprintf(line, " %-+4.4f  %-+4.4f\n", city.lon/deg_to_rad, city.lat/deg_to_rad); 
    outfile << line; 
  }
  outfile.close(); 

  cout << "done." << endl; 
}

//__________________________________________________________________________________________
//returns the total length of a fully-connected path, in radians
double total_length(const vector<CityCoord>& cities, const vector<int>& order) {

  if (cities.size() != order.size()) { throw invalid_argument("order-vector and cities do not match"); }

  double length = 0.; 
  for (size_t i=0; i<cities.size()-1; i++) {

    length += CityDistance( cities[order[i]], cities[order[i+1]] );
  }
  length += CityDistance( cities[order.back()], cities[order.front()] );

  return length; 
}

int main(int argc, char *argv[]){

  int n_test_cities = 200; 
  vector<CityCoord> cities;
  
  //Arg descriptions: 
  // -f   path to list of 'real' cities to use. fake cities will be generated if none are provided
  // -g   if using fake cities, this is the generation technique for the fake cities ('equatorial' or 'random')
  // -n   number of fake cities to generate
  // -d   draw annealing history in plot with this plot name.
  // -o   output final arrangement of cities (arg is file path to put it under)
  NiceOptionParser options(argc, argv, "f:g:n:d:o:");

  bool use_fake_cities = !options.Is_option_set('f'); 

  if (!use_fake_cities) { 

    //get list of real cities 
    char* file_path = options.Arg_or_defalut_value('f', "").data(); 

    printf("Using real city data from file: %s\n", file_path); 

    cities = GetData(file_path); 
      
    printf("Read %zi cities from data file\n", cities.size());  

  } else {

    //generate 'fake' cities randomly
    string gen_type = options.Arg_or_defalut_value('g', "equatorial");

    n_test_cities = atoi(options.Arg_or_defalut_value('n', "200").data()); 
    
    cities = MakeFakeCities( n_test_cities, gen_type );

    printf("Using %i randomly-generated cities with generation rule '%s'\n", n_test_cities, gen_type.data() ); 
  }



  //printf("Longitude  Latitude\n");
  //for (const auto& city : cities) printf("%lf %lf\n",	city.lon, city.lat); 

  //make a random index generator
  random_device rd; 
  mt19937 randgen(rd()); 

  auto rand_index = std::bind( uniform_int_distribution<int>(0, cities.size()-1), std::ref(randgen) ); 

  auto rand_uniform = std::bind( uniform_real_distribution<double>(0., 1.), std::ref(randgen) ); 

  const int n_cities = cities.size(); 

  //create an ordering of cities. 
  // this function picks a random starting city, and goes through each city 1-by-1, picking the closest (remaining)
  // city each time to pair up with. 
  // 
  // obviously, this isn't acceptable as a final answer, but the idea is that its a decent first guess. 
  //
  vector<int> order = SortCities(cities); 

  double length_greedy = total_length(cities, order);
  printf("Greedy-search length: %.2e km  (%.2f circumferences)\n", length_greedy * earth_radius_km, length_greedy/(2.*pi)); 

  //initialize the annealer object 
  Annealer annealer(cities); 
  
  //shuffle the list of cities ('melting')
  annealer.Shuffle(); 
  
  //Structure of Annealing Schedule: 
  // 1. there are a fixed number of 'steps' (n_steps). 
  // 2. For each step, a fixed number of 'swaps' are performed (n_swaps_per_step).
  //    a single 'swap' is a proposed update in which the order of two cities is exchanged. 
  //    the temperature (T) is updated at every step, so it is constant for each swap in the step. 
  // 3. the total length at the end of each step is recorded and plotted. 
  //
  
  //this is the number of times that T is updated 
  const long int n_steps = 250; 

  //starting & ending T-values (they decrease exponentially)
  const double T_start = 2.5;
  const double T_end   = 1.25e-3;  
  //the ratio by which beta increases for each 'T-update step'. 
  const double inflation_per_step = pow(T_end/T_start, 1./((double)n_steps-1)); 

  //I'm assuming that the strength of the annealing algorithm goes as the square of 'n_cities', as 
  const long int n_swaps_per_step = ((long int)n_cities*n_cities) * 200; 

  //these keep track of the 'Temp' and the current total length, so we can print the annealing schedule 
  vector<double> pts_length, pts_temp; 

  //the starting length from the fully 'melted' city list 
  const double starting_length = annealer.Get_total_length() * earth_radius_km; 


  //print out annealing parameter information. 
  printf(
    "\n"
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    "Annealing parameters:\n"
    " T-start:          %.3e\n"
    " T-end:            %.3e\n"
    " n_steps:          %li\n"
    " n_swaps_per_step: %li\n"
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
    T_start, 
    T_end,
    n_steps, 
    n_swaps_per_step
  ); cout << flush; 

  double T = T_start; 

  //use this to time the loop execution
  auto time_start = chrono::steady_clock::now(); 

  for (long int i=0; i<n_steps; i++) {

    //take the step, and measure the time
    auto tstep_start = chrono::steady_clock::now(); 

    //perform the annealing swaps
    double fraction_accepted = annealer.Swap_update(T, n_swaps_per_step);    
    
    //measure the new, total length
    double length = annealer.Get_total_length(); 
    
    //stop the timer 
    auto tstep_end = chrono::steady_clock::now(); 

    chrono::duration<double, std::milli> duration_step = tstep_end - tstep_start; 

    //record the new length (and temperature)
    pts_length.push_back( length * earth_radius_km ); 
    pts_temp  .push_back( T * earth_radius_km ); //our temperature is expressed in radians, so let's convert it to KM 
    printf("\rstep %4li (%3.1f %)  T: %+.2e   length %.2e km  (%.2f circumferences) Fraction of swaps accepted: %.7f  time: %.1f ms", 
      i, 
      100.*((double)i)/((double)n_steps), 
      T * earth_radius_km, 
      length * earth_radius_km, 
      length/(2.*pi), 
      fraction_accepted, 
      duration_step.count()
    ); cout << flush; 
    
    //update the temperature 
    T *= inflation_per_step; 
  }

  auto time_end = chrono::steady_clock::now(); 

  
  //time elapsed, in seconds
  chrono::duration<double, std::milli> execution_duration = time_end - time_start; 
  const double time_elapsed_us = execution_duration.count(); 


  //final length after annealing is done
  const double ending_length = annealer.Get_total_length() * earth_radius_km; 

  printf(
    "\n"
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    "Iterations done. \n"
    " - time elapsed:     %.3f s (%.1f ms / step)\n"
    "\n"
    " - starting length:  %.1f x 10^3 km  (%.4f x earth circumference)\n"
    " - ending length:    %.1f x 10^3 km  (%.4f x earth circumference)\n"
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", 
    time_elapsed_us/1.e3, time_elapsed_us/((double)n_steps),

    starting_length/1.e3, starting_length/(2.*pi*earth_radius_km), 
    ending_length/1.e3, ending_length/(2.*pi*earth_radius_km) 
  );


  //save the new city-path, if that option is specified. 
  if (options.Is_option_set('o')) {
    
    string path_out = options.Arg_or_defalut_value('o', "cities_out.dat"); 
    save_city_csv(path_out, cities, annealer.Get_city_order()); 
  }

  auto c = new TCanvas; 

  gPad->SetLogx(1); 
  auto g = new TGraph(pts_temp.size(), pts_temp.data(), pts_length.data()); 
  gPad->SetLeftMargin(0.15); 
  g->SetTitle("Annealing Schedule;Temperature (units: km);total length (units: km)"); 
  g->GetYaxis()->SetRangeUser(0., vector_max(pts_length)*1.12 ); 
  g->Draw(); 

  string plot_name = options.Arg_or_defalut_value('d', "annealing_test.png"); 
  c->SaveAs(plot_name.data()); 

  return 0;
}


// example code to read in a data file of city lat,long coordinates

#include <cstdio>
#include <cstdlib>
#include <vector>
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
#include <stdexcept> 
#include <fstream> 

using namespace std; 

namespace {
  constexpr double pi = 3.14159265359; 

  
}
//get the maximum element from a std::vector<T>
template<typename T> T vector_max(const std::vector<T>& V, T max=(T)0.) {
  for (const T& x : V) if (x > max) max = x; 
  return max; 
}

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
//returns the total length of a fully-connected path, in radians
void save_city_csv(const string& path_out, const vector<CityCoord>& cities, const vector<int>& order)
{
  cout << "\nMaking output file '" << path_out << "'..." << flush; 

  //create the output file 
  ofstream outfile(path_out.data(), ios::trunc); 

  outfile << "# longitude latitude\n"; 
  for (int index : order) {
    const auto& city = cities[index];
    
    char line[500]; sprintf(line, " %-+4.4f  %-+4.4f\n", city.lon, city.lat); 
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

    char* file_path = options.Arg_or_defalut_value('f', "").data(); 

    printf("Using real city data from file: %s\n", file_path); 

    cities = GetData(file_path); 

  } else {
    
    string gen_type = options.Arg_or_defalut_value('g', "equatorial");

    n_test_cities = atoi(options.Arg_or_defalut_value('n', "200").data()); 
    
    cities = MakeFakeCities( n_test_cities, gen_type );

    printf("Using %i randomly-generated cities with generation rule '%s'\n", n_test_cities, gen_type.data() ); 
  }

  printf("Read %zi cities from data file\n", cities.size()); 


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
#if 0 
  //_________________________________________________________________________________________________________________
  auto annealing_update = [&rand_index,&rand_uniform,&cities,n_cities](double Beta, vector<int>& order)
  {
    //check the 'energy change' (change in total angle) between the current state and the proposed update. 

    //get two (distinct!) random indices
    int ind1 = rand_index(); 
    int ind2; do { ind2 = rand_index(); } while (ind1==ind2);

    
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
        CityDistance( 
          cities[i_lo-1>=0 ? order[i_lo-1] : order.back()], 
          cities[order[i_hi]] 
        ) + 
        CityDistance( 
          cities[i_hi+1<n_cities ? order[i_hi+1] : order.front()], 
          cities[order[i_lo]] 
        ); 

      //now, compute what these distances are in the original (current) ordering. 
      d_energy -= 
        CityDistance( 
          cities[i_lo-1>=0 ? order[i_lo-1] : order.back()], 
          cities[order[i_lo]] 
        ) + 
        CityDistance( 
          cities[i_hi+1<n_cities ? order[i_hi+1] : order.front()], 
          cities[order[i_hi]] 
        ); 

    } else { //this is for cities that are NOT currently adjacent. 
        
      //this only works for NON-ADJACENT cities (adjaceny defined by the current ordering of cities, NOT by physical proximity)
      auto check_swapped_distances = [&](int index_old, int index_new)
      {
        return
          CityDistance( 
            cities[index_old-1>=0 ? order[index_old-1] : order.back()], 
            cities[order[index_new]] 
          ) +  
          CityDistance( 
            cities[order[index_new]], 
            cities[index_old+1<n_cities ? order[index_old+1] : order.front()] 
          ); 
      };

      d_energy = 
        check_swapped_distances(ind1, ind2) + 
        check_swapped_distances(ind2, ind1); 

      d_energy -= 
        check_swapped_distances(ind1, ind1) + 
        check_swapped_distances(ind2, ind2); 
    }

    double prob_increase = exp( -Beta*d_energy ); 

    if (d_energy < 0. || prob_increase > rand_uniform()) {

      //then, accept the new update. swap the orders of the cities
      int index_new_1 = order[ind2]; 
      int index_new_2 = order[ind1]; 

      order[ind1] = index_new_1;  
      order[ind2] = index_new_2; 

      return true; //new update accepted
    }
    return false; //new update rejected. 
  };
  //_________________________________________________________________________________________________________________
#endif 

  //starting beta value
  double beta = 2.; 

  //initialize the annealer object 
  Annealer annealer(cities); 
  
  const long int n_steps = 8e8; 

  const long int report_step = n_steps / 500; 

  //the ratio by which beta increases for each step
  const double inflation_per_step = pow(50./beta, 1./((double)n_steps)); 

  vector<double> len, step_over_N; 

  for (long int i=0; i<n_steps; i++) {
    
    annealer.Swap_update(beta); 

    beta = 1. + 25.*( pow( ((double)n_steps-i)/((double)n_steps), -0.5 ) - 1 ); 

    if (i % report_step == 0) {

      double length = annealer.Get_total_length(); 

      len.push_back( length / (2.*pi) ); 
      step_over_N.push_back(((double)i)/((double)n_cities)); 
      printf("\rstep %8li   beta %+.2e    length (in earth circumferences:) %4.6f", i, beta, length/(2.*pi)); cout << flush; 
    }
  }


  //save the new city-path, if that option is specified. 
  if (options.Is_option_set('o')) {
    
    string path_out = options.Arg_or_defalut_value('o', "cities_out.dat"); 
    save_city_csv(path_out, cities, order); 
  }

  auto c = new TCanvas; 

  //gPad->SetLogx(1); 
  auto g = new TGraph(len.size(), step_over_N.data(), len.data()); 
  gPad->SetLeftMargin(0.15); 
  g->SetTitle("Total length vs step;step / N. cities;total length (units: Earth Circumference)"); 
  g->GetYaxis()->SetRangeUser(0., vector_max(len)*1.12 ); 
  g->Draw(); 

  string plot_name = options.Arg_or_defalut_value('d', "annealing_test.png"); 
  c->SaveAs(plot_name.data()); 

  return 0;
}


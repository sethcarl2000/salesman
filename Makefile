default: all


all: iterate simulated_annealing benchmark

iterate: iterate.cpp simulated_annealing
	g++ -Wall iterate.cpp simulated_annealing.o -o iterate

iterate: benchmark.cpp simulated_annealing
	g++ -Wall benchmark.cpp simulated_annealing.o -o benchmark

simulated_annealing: simulated_annealing.hpp CityCoord.hpp simulated_annealing.cpp CityCoord.hpp
	g++ -Wall simulated_annealing.cpp -c simulated_annealing.o


clean:
	rm -f iterate *~ *png *.o 

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  += $(ROOTCFLAGS) -Wall -fPIC
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) # ,-rpath #,$(ODELIB) -L$(ODELIB) -lode
GSLFLAGS   = -lgsl -lgslcblas
GXX	   = g++ $(CXXFLAGS)


default: all

all: iterate simulated_annealing benchmark libNiceOptionParser.so

# compile the library into an object file 
libNiceOptionParser.o: NiceOptionParser.hpp NiceOptionParser.cpp
	$(GXX) -c NiceOptionParser.cpp -o libNiceOptionParser.o

# Create the shared library
libNiceOptionParser.so: libNiceOptionParser.o
	$(GXX) -shared -o libNiceOptionParser.so libNiceOptionParser.o


iterate: iterate.cpp simulated_annealing libNiceOptionParser.so
	$(GXX) iterate.cpp simulated_annealing.o -o iterate -L. -lNiceOptionParser -Wl,-rpath,'$$ORIGIN' $(ROOTFLAGS) 

simulated_annealing: simulated_annealing.hpp CityCoord.hpp simulated_annealing.cpp CityCoord.hpp libNiceOptionParser.so
	$(GXX) -Wall simulated_annealing.cpp -c simulated_annealing.o $(ROOTFLAGS) -L. -lNiceOptionParser


clean:
	rm -f iterate *~ *png *.o *.so

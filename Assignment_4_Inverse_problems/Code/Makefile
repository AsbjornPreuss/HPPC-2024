CXX := clang++-12
CXX := g++ 

#######################################################
# Optimization flags are chosen as the last definition.
# Comment out using "#" at the begining of the line or rearrange according to your needs.
#
# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT=-O3 -ffast-math

# Add thread checking to code
#OPT=-O1 -fsanitize=thread -fPIE -pie

# Faster compilation time
#OPT=-O1

#Profiling
OPT=-O1 -pg

CXXFLAGS := $(OPT) -Wall -march=native -g -std=c++14

default: seq mp

seq: seismogram_seq.cpp
	$(CXX) seismogram_seq.cpp $(CXXFLAGS) -o seq

mp: seismogram_omp_l3.cpp
	$(CXX) seismogram_omp_l3.cpp $(CXXFLAGS) -fopenmp -o mp_3_4

clean:
	rm -fr seq mp

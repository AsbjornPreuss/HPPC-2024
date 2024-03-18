LIB = -L/opt/openmpi-4.1.6/lib -Wl,-rpath -Wl,/opt/openmpi-4.1.6/lib -Wl,--enable-new-dtags -lmpi
INC = -I/opt/openmpi-4.1.6/include

CXX := mpic++

FLAGS = -O1 -DNDEBUG -g -Wall

.PHONY: clean all

all: seq_hsb parallel_hsb

seq_hsb: seq_hsb.cpp
	$(CXX) $(FLAGS) $(LIB) $(INC) -o seq_hsb seq_hsb.cpp

parallel_hsb: parallel_hsb.cpp
	$(CXX) $(FLAGS) $(LIB) $(INC) -o parallel_hsb parallel_hsb.cpp

clean:
	rm -f seq_hsb parallel_hsb *.o

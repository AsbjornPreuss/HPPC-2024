FLAGS = -O3 -DNDEBUG -g -Wall

.PHONY: clean all

all: seq_hsb parallel_hsb

seq_hsb: seq_hsb.cpp
	g++ $(FLAGS) $(LIB) $(INC) -o seq_hsb seq_hsb.cpp

parallel_hsb: parallel_hsb.cpp
	g++ $(FLAGS) $(LIB) $(INC) -o parallel_hsb parallel_hsb.cpp

clean:
	rm -f seq_hsb parallel_hsb *.o

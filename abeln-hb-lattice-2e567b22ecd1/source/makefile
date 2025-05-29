PROGRAM=hblattice
PROGNAT=grancanNat
#CC=mpiCC
#CC=mpic++
CC=mpicxx
CFLAGS=-O3 -Wall 
SOURCES=ptemp.cpp Chain.cpp Lattice.cpp Stats.cpp EnergyMap.cpp Native.cpp AA.cpp MolBox.cpp GranCan.cpp Cluster.cpp Pos.cpp Design.cpp App.cpp AADistr.cpp MonteCarlo.cpp Moves.cpp
OBJECTS=$(SOURCES:.cpp=.o)
HEADERS=ptemp.h Lattice.hpp Chain.hpp Residue.hpp Stats.hpp EnergyMap.hpp Native.hpp AA.hpp MolBox.hpp GranCan.hpp Cluster.hpp Pos.hpp StatsMethods.hpp Design.hpp App.hpp AADistr.hpp MonteCarlo.hpp Moves.hpp

all: $(PROGRAM) 

$(PROGRAM): $(OBJECTS) 
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) 


#
# Object targets: rules that define objects, their dependencies, and
# a list of commands for compilation.
#

#$(OBJECTS): $(SOURCES) $(HEADERS)
#	$(CC) $(CFLAGS) -c -o $(OBJECTS) $(HEADERS)


%.o : %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $<


# Clean target: "make clean" to remove unwanted objects and executables.
#

clean:
	rm $(OBJECTS) $(PROGRAM) 

#
# Run target: "make -f Makefile.Linux run" to execute the application
#             You will need to add $(VARIABLE_NAME) for any command line parameters 
#             that you defined earlier in this file.
# 

run:
	./$(PROGRAM) 


tidy: 
	rm *~

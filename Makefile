# To silence PMIX errors: export PMIX_MCA_gds=hash

CXX = mpicxx
CXXFLAGS = -std=c++14 -Wall -O3


SOURCES = src/definitions.cpp \
          src/lattice.cpp \
          src/ising.cpp \
          src/mcrg.cpp \
          src/main.cpp
          
OBJECTS = build/definitions.o \
          build/lattice.o \
          build/ising.o \
          build/mcrg.o \
          build/main.o
          
          
all: run


clean:
	rm -f run ${OBJECTS} slurm*

run: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o run
	
build/main.o: src/main.cpp
	${CXX} ${CXXFLAGS} -c src/main.cpp -o build/main.o

build/definitions.o: src/definitions.cpp
	${CXX} ${CXXFLAGS} -c src/definitions.cpp -o build/definitions.o

build/lattice.o: src/lattice.cpp
	${CXX} ${CXXFLAGS} -c src/lattice.cpp -o build/lattice.o

build/ising.o: src/ising.cpp
	${CXX} ${CXXFLAGS} -c src/ising.cpp -o build/ising.o
	
build/mcrg.o: src/mcrg.cpp
	${CXX} ${CXXFLAGS} -c src/mcrg.cpp -o build/mcrg.o

# To silence PMIX errors: export PMIX_MCA_gds=hash

CXX = mpicxx
CXXFLAGS = -std=c++14 -Wall -O3


SOURCES = src/definitions.cpp \
          src/lattice.cpp \
          src/ising.cpp \
          src/mcrg.cpp \
          src/rgnn.cpp \
          
OBJECTS = build/definitions.o \
          build/lattice.o \
          build/ising.o \
          build/mcrg.o \
          build/rgnn.o
          

all: run train


clean:
	rm -f run train build/* jobs/* slurm*

run: ${OBJECTS} build/main.o
	${CXX} ${CXXFLAGS} ${OBJECTS} build/main.o -o run

train: ${OBJECTS} build/train.o
	${CXX} ${CXXFLAGS} ${OBJECTS} build/train.o -o train
	
build/main.o: src/main.cpp
	${CXX} ${CXXFLAGS} -c src/main.cpp -o build/main.o

build/train.o: src/train.cpp
	${CXX} ${CXXFLAGS} -c src/train.cpp -o build/train.o

build/definitions.o: src/definitions.cpp
	${CXX} ${CXXFLAGS} -c src/definitions.cpp -o build/definitions.o

build/lattice.o: src/lattice.cpp
	${CXX} ${CXXFLAGS} -c src/lattice.cpp -o build/lattice.o

build/ising.o: src/ising.cpp
	${CXX} ${CXXFLAGS} -c src/ising.cpp -o build/ising.o
	
build/mcrg.o: src/mcrg.cpp
	${CXX} ${CXXFLAGS} -c src/mcrg.cpp -o build/mcrg.o
	
build/rgnn.o: src/rgnn.cpp
	${CXX} ${CXXFLAGS} -c src/rgnn.cpp -o build/rgnn.o

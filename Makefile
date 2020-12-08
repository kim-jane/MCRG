# To silence PMIX errors: export PMIX_MCA_gds=hash

CXX = mpicxx
CXXFLAGS = -std=c++14 -Wall -O3


SOURCES = src/definitions.cpp \
          src/ising.cpp \
          src/mcrg.cpp \
          src/test_ising.cpp \
          src/test_mcrg.cpp \
          src/main.cpp
          
OBJECTS = build/definitions.o \
          build/ising.o \
          build/mcrg.o \
          build/test_ising.o \
          build/test_mcrg.o \
          build/main.o
          
          
all: run


clean:
	rm -f run ${OBJECTS}

run: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o run
	
build/main.o: src/main.cpp
	${CXX} ${CXXFLAGS} -c src/main.cpp -o build/main.o

build/definitions.o: src/definitions.cpp
	${CXX} ${CXXFLAGS} -c src/definitions.cpp -o build/definitions.o

build/ising.o: src/ising.cpp
	${CXX} ${CXXFLAGS} -c src/ising.cpp -o build/ising.o
	
build/mcrg.o: src/mcrg.cpp
	${CXX} ${CXXFLAGS} -c src/mcrg.cpp -o build/mcrg.o
	
build/test_ising.o: src/test_ising.cpp
	${CXX} ${CXXFLAGS} -c src/test_ising.cpp -o build/test_ising.o

build/test_mcrg.o: src/test_mcrg.cpp
	${CXX} ${CXXFLAGS} -c src/test_mcrg.cpp -o build/test_mcrg.o



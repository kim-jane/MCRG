CXX = g++
CXXFLAGS = -std=c++14 -Wall -O3


SOURCES = definitions.cpp \
          ising.cpp \
          mcrg.cpp \
          main.cpp
          
OBJECTS = definitions.o \
          ising.o \
          mcrg.o \
          main.o
          
          
all: run


clean:
	rm -f run ${OBJECTS}

run: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o run

definitions.o: definitions.cpp
	${CXX} ${CXXFLAGS} -c definitions.cpp -o definitions.o

ising.o: ising.cpp
	${CXX} ${CXXFLAGS} -c ising.cpp -o ising.o
	
mcrg.o: mcrg.cpp
	${CXX} ${CXXFLAGS} -c mcrg.cpp -o mcrg.o

main.o: main.cpp
	${CXX} ${CXXFLAGS} -c main.cpp -o main.o

CXX = g++
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
	

definitions.o: definitions.cpp
	${CXX} ${CXXFLAGS} -c definitions.cpp -o definitions.o

ising.o: ising.cpp
	${CXX} ${CXXFLAGS} -c ising.cpp -o ising.o
	
mcrg.o: mcrg.cpp
	${CXX} ${CXXFLAGS} -c mcrg.cpp -o mcrg.o
	
test_ising.o: test_ising.cpp
	${CXX} ${CXXFLAGS} -c test_ising.cpp -o test_ising.o

test_mcrg.o: test_mcrg.cpp
	${CXX} ${CXXFLAGS} -c test_mcrg.cpp -o test_mcrg.o

main.o: main.cpp
	${CXX} ${CXXFLAGS} -c main.cpp -o main.o

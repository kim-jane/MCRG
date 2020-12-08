#include "mcrg.hpp"
#include "test_ising.hpp"


int main(){
    
    MPI_Init(NULL, NULL);

    TestIsing2D test;
    test.equilibrate_sweepK1(1E6, 64, 0.4, 0.5, 0.01);

    MPI_Finalize();
    return 0;
}

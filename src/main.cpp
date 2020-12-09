#include "mcrg.hpp"
#include "test_ising.hpp"


int main(){
    
    MPI_Init(NULL, NULL);

    TestIsing2D test;
    test.equilibrate_sweepK1(1E5, 16, -0.50, -0.40, 0.01);

    
    MPI_Finalize();
    return 0;
}

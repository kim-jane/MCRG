#include "mcrg.hpp"
#include "test_ising.hpp"


int main(){
    
    MPI_Init(NULL, NULL);

    TestIsing2D test;
    test.equilibrate_sweepK1(1E6, 16, -1.0, -0.2, 0.1);
    test.equilibrate_sweepK1(1E6, 32, -1.0, -0.2, 0.1);
    test.equilibrate_sweepK1(1E6, 64, -1.0, -0.2, 0.1);
    
    MPI_Finalize();
    return 0;
}

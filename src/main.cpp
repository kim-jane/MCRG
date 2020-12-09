#include "mcrg.hpp"
#include "test_ising.hpp"


int main(){
    
    MPI_Init(NULL, NULL);

    TestIsing2D test;
    test.equilibrate_sweepK1(2E5, 4, -1.0, -0.2, 0.1);
    test.equilibrate_sweepK1(2E5, 128, -1.0, -0.2, 0.1);
    
    /*
    int b = 2;
    vec2D Kc;
    Kc(0) = -log(1+sqrt(2))/2;
    Kc(1) = 0.0;
    
    MonteCarloRenormalizationGroup MCRG(b);
    double nu = MCRG.calc_critical_exponent(1E6, 64, Kc);
    
    */
    
    MPI_Finalize();
    return 0;
}

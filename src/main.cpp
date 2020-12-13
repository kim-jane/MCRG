#include "mcrg.hpp"
#include "test_ising.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    vec2D Kc;
    Kc(0) = -0.45;
    Kc(1) = 0.0;

    /*
    TestIsing2D test;
    test.check_block_spin_transformation(1E7, b, 128, Kc);
    test.equilibrate_sweepT(1E4, 64, 1.0, 4.0, 0.1);
    */

    
    MonteCarloRenormalizationGroup MCRG(b);
    
    Kc = MCRG.locate_critical_point(100, 1E4, 1E6, 64, Kc);
    MCRG.calc_critical_exponent(10, 1E5, 1E6, 64, Kc);
    
    MPI_Finalize();
    return 0;
}

#include "mcrg.hpp"
#include "test_ising.hpp"


int main(){
    
    MPI_Init(NULL, NULL);

    int b = 2;
    vec2D Kc;
    Kc(0) = -log(1+sqrt(2))/2;
    Kc(1) = 0.0;
    
    MonteCarloRenormalizationGroup MCRG(b);
    //double nu = MCRG.calc_critical_exponent(1E6, 32, Kc);
    
    vec2D Kc_approx = MCRG.approx_critical_point(1E6, 32, Kc);
    
    MPI_Finalize();
    return 0;
}

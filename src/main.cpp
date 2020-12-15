#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    
    int b = 2;
    vec2D Kc;
    Kc(0) = -0.44;
    Kc(1) = 0.0;
     

    MonteCarloRenormalizationGroup MCRG(b);
    
    Kc = MCRG.locate_critical_point(10, 1E4, 1E6, 64, Kc);
    //MCRG.calc_critical_exponent(1, 1E4, 1E6, 32, Kc);
    
    MPI_Finalize();
    return 0;
}

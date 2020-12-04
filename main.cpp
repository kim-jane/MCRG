#include "mcrg.hpp"


int main(){

    int b = 3;           // scaling factor
    int N0 = 81;         // initial lattice size
    int n_samples = 1E6; // MC samples
    
    vec2D K0;            // initial couplings
    K0[0] = -1.0;        // nearest neighbor
    K0[1] = 0.0;         // next nearest neighbor
    
    // known critical temperature
    double Tc = -2.0*K0[0]/log(1.0+sqrt(2.0));
    //double Tc = 1.0;
    
    MonteCarloRenormalizationGroup MCRG(b, 1.0);
    MCRG.locate_critical_point(n_samples, N0, K0*Tc);
    //MCRG.calc_critical_exponent();
    
    return 0;
}

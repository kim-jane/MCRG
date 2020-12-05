#include "mcrg.hpp"


int main(){

    int b = 3;           // scaling factor
    int N0 = 81;         // initial lattice size
    int n_samples = 1E6; // MC samples
    
    vec2D K0;            // initial couplings
    K0[0] = -log(1+sqrt(2))/2;        // nearest neighbor
    K0[1] = 0.0;         // next nearest neighbor
    
    MonteCarloRenormalizationGroup MCRG(b);
    MCRG.run(n_samples, N0, K0);

    
    return 0;
}

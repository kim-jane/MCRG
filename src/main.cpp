#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);

    int b = 2;               // scaling factor
    int L0 = 128;
    //int N0 = 64;             // initial lattice size
    
    int n_samples = 1E6;     // MC samples
    int verbose = 1;
    
    vec2D K0;                // initial couplings
    K0[0] = -0.44;           // nearest neighbor
    K0[1] = 0.0;             // next nearest neighbor
    
    //double nu;
    vec2D Kc;

    MonteCarloRenormalizationGroup MCRG(b, verbose);
    Kc = MCRG.locate_critical_point(10, n_samples, L0, K0);
    //nu = MCRG.calc_critical_exponent(10*n_samples, 2*N0, Kc);
    

    MPI_Finalize();
    return 0;
}

#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 16;
    int n_samples_eq = 1E4;
    int n_samples = 1E7;
    //double K = -0.44068277273;
    double Kc = -log(1+sqrt(2))/2;

    MonteCarloRenormalizationGroup MCRG(b);

    //MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, K);
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc);
    
    MPI_Finalize();
    return 0;
}

/*
Average K for different lattice sizes:
16 -0.440414806102000
32 -0.440619259540000
64 -0.440675235450000
128 -0.440682772730000
*/

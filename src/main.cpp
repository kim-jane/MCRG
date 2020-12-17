#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 32;
    int n_samples_eq = 1E4;
    int n_samples = 1E6;
    double K = -0.440619;

    MonteCarloRenormalizationGroup MCRG(b);

    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, K);
    
    MPI_Finalize();
    return 0;
}

#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 16;
    int n_samples_eq = 1E4;
    int n_samples = 1E6;
    int n_iterations = 80;
    double K = -0.440415;

    MonteCarloRenormalizationGroup MCRG(b);

    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, K);
    
    MPI_Finalize();
    return 0;
}

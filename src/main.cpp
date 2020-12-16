#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 32;
    int n_samples_eq = 1E4;
    int n_samples = 1E6;
    int n_iterations = 10;
    double K = -0.44;

    MonteCarloRenormalizationGroup MCRG(b);
    MCRG.locate_critical_point(n_iterations, n_samples_eq, n_samples, N, K);
    MCRG.calc_critical_exponent(1E4, 1E6, 64, K);
    
    MPI_Finalize();
    return 0;
}

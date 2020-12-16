#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 128;
    int n_samples_eq = 1E4;
    int n_samples = 1E6;
    int n_iterations = 20;
    double K0 = -0.44;
    double Kc = -log(1+sqrt(2))/2;

    MonteCarloRenormalizationGroup MCRG(b);
    
    // calc thermal exponent at critical coupling for N->inf
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc);
    
    // locate critical point for finite N
    double Kc_N = MCRG.locate_critical_point(n_iterations, n_samples_eq, n_samples, N, K0);
    
    // calc thermal exponent at critical point
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc_N);
    
    
    MPI_Finalize();
    return 0;
}

#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    
    int n_samples_eq = 1E4;
    int n_samples = 1E6;
    int n_iterations = 100;
    double K0 = -0.44;
    double Kc = -log(1+sqrt(2))/2;
    double Kc_N;

    MonteCarloRenormalizationGroup MCRG(b);

    int N = 32;
    
    // calc thermal exponent at critical coupling for N->inf
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc);
    
    // locate critical point for finite N
    Kc_N = MCRG.locate_critical_point(n_iterations, n_samples_eq, n_samples, N, K0);
    
    // calc thermal exponent at critical point
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc_N);
    
    N = 64;
    
    // calc thermal exponent at critical coupling for N->inf
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc);
    
    // locate critical point for finite N
    Kc_N = MCRG.locate_critical_point(n_iterations, n_samples_eq, n_samples, N, K0);
    
    // calc thermal exponent at critical point
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc_N);
    
    N = 128;
    
    // calc thermal exponent at critical coupling for N->inf
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc);
    
    // locate critical point for finite N
    Kc_N = MCRG.locate_critical_point(n_iterations, n_samples_eq, n_samples, N, K0);
    
    // calc thermal exponent at critical point
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc_N);
    
    
    MPI_Finalize();
    return 0;
}

#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 64;
    int n_samples_eq = 1E4;
    int n_samples = 1E6;
    int n_iterations = 1000;
    double K0 = -0.44;
    double Kc = -log(1+sqrt(2))/2;

    MonteCarloRenormalizationGroup MCRG(b);
    
    double Kc_N = MCRG.locate_critical_point(n_iterations, n_samples_eq, n_samples, N, K0);
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc_N);
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc);
    
    MPI_Finalize();
    return 0;
}

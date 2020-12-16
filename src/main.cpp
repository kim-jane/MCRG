#include "mcrg.hpp"


int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 64;
    int n_samples_eq = 1E4;
    int n_samples = 1E7;
    int n_iterations = 10;
    double K0 = -0.40;
    double Kc = -log(1+sqrt(2))/2;

    MonteCarloRenormalizationGroup MCRG(b);
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc);
    
    double Kc_N = MCRG.locate_critical_point(n_iterations, n_samples_eq, n_samples, N, K0);
    MCRG.calc_critical_exponent(n_samples_eq, n_samples, N, Kc_N);
    
    
    MPI_Finalize();
    return 0;
}

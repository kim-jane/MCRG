#include "rgnn.hpp"

int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 32;
    int n_cycles = 1E4;
    int n_samples = 1E5;
    int n_samples_eq = 1E4;
    double Tc = 2/log(1+sqrt(2));
    double h = 0.0001;
    double eta = 0.001;
    
    RenormalizationGroupNeuralNetwork RGNN(b);
    
    RGNN.train(N, n_cycles, n_samples, n_samples_eq, Tc, h, eta);
    
    MPI_Finalize();
    return 0;
}

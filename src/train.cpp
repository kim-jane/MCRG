#include "rgnn.hpp"

int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 64;
    int n_cycles = 1E3;
    int n_samples = 1E6;
    int n_samples_eq = 1E4;
    double Tc = 1.0;
    double h = 0.001;
    double eta = 0.001;
    
    RenormalizationGroupNeuralNetwork RGNN(b);
    
    RGNN.train(N, n_cycles, n_samples, n_samples_eq, Tc, h, eta);
    
    MPI_Finalize();
    return 0;
}

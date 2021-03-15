#include "rgnn.hpp"

int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = 128;
    int n_cycles = 1E4;
    int n_samples = 1E6;
    int n_samples_eq = 1E3;
    double K = -log(1+sqrt(2))/2;
    
    RenormalizationGroupNeuralNetwork RGNN(b);
    
    RGNN.train(N, n_cycles, n_samples, n_samples_eq, K);
    
    MPI_Finalize();
    return 0;
}

#include "rgnn.hpp"

int main(int argc, char* argv[]){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int n_cycles = 1E3;
    int n_samples = 1E4;
    int n_samples_eq = 1E4;
    double Kc = -log(1+sqrt(2))/2;
    double h = 0.0001;
    double eta = 0.001;
    
    
    RenormalizationGroupNeuralNetwork RGNN(b);
    
    RGNN.train_scalar_output(8, 10*n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    
    for(int N = 16; N <= 64; N *= b){
        RGNN.train_scalar_output(N, n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    }
    
    for(int N = 8; N <= 64; N *= b){
        RGNN.test_scalar_output(N, 10*n_samples, n_samples_eq, Kc, 0.1);
    }
    
    MPI_Finalize();
    return 0;
}

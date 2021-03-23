#include "rgnn.hpp"

int main(int argc, char* argv[]){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int n_cycles = 1E3;
    int n_samples = 1E4;
    int n_samples_eq = 1E3;
    double Kc = -log(1+sqrt(2))/2;
    double h = 0.0001;
    double eta = 0.001;
    
    // initialize network with filter size b
    RenormalizationGroupNeuralNetwork RGNN(b);
    
    
    // SUPERVISED LEARNING

    // train output to predict temperature at Kc
    // start with smallest input size for many cycles
    RGNN.train_temperature(8, 10*n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    
    // train with larger input sizes for fewer cycles
    for(int N = 16; N <= 64; N *= b){
        RGNN.train_temperature(N, n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    }
    
    // test network in range around Kc
    for(int N = 8; N <= 64; N *= b){
        RGNN.test_temperature(N, 100*n_samples, 100*n_samples_eq, Kc, 0.4);
    }
    
    
    // reinitialize network
    RGNN.initialize();

    // UNSUPERVISED LEARNING
    
    // train output to be scale invariant at Kc
    // start with smallest input size for many cycles
    RGNN.train_scalar_output(8, 10*n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    
    // train with larger input sizes for fewer cycles
    for(int N = 16; N <= 64; N *= b){
        RGNN.train_scalar_output(N, n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    }
    
    // test network in range around Kc
    for(int N = 8; N <= 64; N *= b){
        RGNN.test_scalar_output(N, 100*n_samples, 100*n_samples_eq, Kc, 0.4);
    }
    
    
    MPI_Finalize();
    return 0;
}

#include "rgnn.hpp"

int main(int argc, char* argv[]){
    
    MPI_Init(NULL, NULL);
    
    int b = atoi(argv[1]);
    int N = atoi(argv[2]);
    int n_cycles = 1E4;
    int n_samples = 1E5;
    int n_samples_eq = 1E4;
    double T = atof(argv[3]);
    double h = 0.0001;
    double eta = 0.001;
    double lambda = 0.01;
    
    RenormalizationGroupNeuralNetwork RGNN(b);
    RGNN.train(N, n_cycles, n_samples, n_samples_eq, T, h, eta, lambda);
    
    MPI_Finalize();
    return 0;
}

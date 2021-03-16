#include "rgnn.hpp"

int main(int argc, char* argv[]){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int N = atoi(argv[1]);
    int n_cycles = 2E3;
    int n_samples = 1E6;
    int n_samples_eq = 1E4;
    double T = atof(argv[2]);
    double h = 0.0001;
    double eta = 0.01;
    
    RenormalizationGroupNeuralNetwork RGNN(b);
      
    RGNN.train(N, n_cycles, n_samples, n_samples_eq, T, h, eta);


    MPI_Finalize();
    return 0;
}

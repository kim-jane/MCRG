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
    
    std::cout << "N = 8" << std::endl;
    RGNN.train_scalar_output(8, 10*n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    
    for(int N = 16; N <= 512; N *= b){
        std::cout << "N = " << N << std::endl;
        RGNN.train_scalar_output(N, n_cycles, n_samples, n_samples_eq, Kc, h, eta);
    }
    
    MPI_Finalize();
    return 0;
}

#include "rgnn.hpp"

int main(){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int n_cycles = 1E3;
    int n_samples = 1E5;
    int n_samples_eq = 1E4;
    double h = 0.0001;
    double eta = 0.001;
    
    // range of lattice sizes
    int Ni = 8;
    int Nf = 128;
    
    // range of temperatures around critical temperature
    double Tc = 2/log(1+sqrt(2));
    double DeltaT = 0.1;
    double dT = 0.01;
    
    RenormalizationGroupNeuralNetwork RGNN(b);
    
    for(int N = Ni; N <= Nf; N *= b){
        for(double T = Tc-DeltaT; T <= Tc+DeltaT; T += dT){
            
            RGNN.initialize_weights();
            RGNN.train(N, n_cycles, n_samples, n_samples_eq, T, h, eta);
        }
    }

    MPI_Finalize();
    return 0;
}

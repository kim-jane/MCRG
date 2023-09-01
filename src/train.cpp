#include "rgnn.hpp"

int main(int argc, char* argv[]){
    
    MPI_Init(NULL, NULL);
    
    int b = 2;
    int n_cycles = 1E4;
    int n_samples = 1E4;
    int n_samples_eq = 1E3;
    //double K0 = -0.3;
    double Tc = 2/log(1+sqrt(2));
    double Kc = -1/Tc;
    double h = 0.0001;
    double eta = 0.001;
    
    int N = 8;
    
    mat W0(2,2);
    W0(0,0) = 0.5;
    W0(0,1) = -0.5;
    W0(1,0) = 0.5;
    W0(1,1) = -0.5;
    
    for(double T = 0.9*Tc; T <= 1.1*Tc; T += 0.2*Tc/5){
        
        double K = -1.0/T;
        
        RenormalizationGroupNeuralNetwork RGNN(b);
        RGNN.set_weights(W0);
        RGNN.train_scalar_output(N, n_cycles, n_samples, n_samples_eq, K, h, eta);
        
        printf("%e %e\n", T, RGNN.final_mse_);
    }
    
    
    /*
    // cost and gradient of cost
    double cost, grad;

    // UNSUPERVISED LEARNING: locating critical nearest neighbor coupling
    int N = 8;
    double K = K0;
    double dK = 0.1;
    for(int cycles = 0; cycles < 10; ++cycles){
        
        // initialize three network with filter size b
        RenormalizationGroupNeuralNetwork RGNN0(b);
        RenormalizationGroupNeuralNetwork RGNN1(b);
        RenormalizationGroupNeuralNetwork RGNN2(b);
        
        // sync weights
        RGNN1.set_weights(RGNN0.W_);
        RGNN2.set_weights(RGNN0.W_);
        
        // train each at different couplings
        RGNN0.train_scalar_output(N, n_cycles, n_samples, n_samples_eq, K, h, eta);
        RGNN1.train_scalar_output(N, n_cycles, n_samples, n_samples_eq, K+dK, h, eta);
        RGNN2.train_scalar_output(N, n_cycles, n_samples, n_samples_eq, K-dK, h, eta);
        
        // estimate gradient of final MSE WRT K
        cost = RGNN0.final_mse_;
        grad = (RGNN1.final_mse_-RGNN2.final_mse_)/(2*dK);
        
        // update coupling
        K -= eta*grad;
        MPI_Bcast(&K, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        printf("C = %lf, grad = %lf, K = %lf\n", cost, grad, K);
    }
    
     */
    
    MPI_Finalize();
    return 0;
}

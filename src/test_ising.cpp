#include "test_ising.hpp"

void TestIsing2D::equilibrate_sweepT(int n_samples_eq,
                                     int N,
                                     double T0,
                                     double Tf,
                                     double dT){
    
    double eps = 1E-5;
    for(double T = T0; T <= Tf+eps; T += dT){
        
        vec2D K;
        K(0) = -1.0/T;
        K(1) = 0.0;
        
        Ising2D* pIsing = new Ising2D(N, K);
        pIsing->equilibrate(n_samples_eq);
    }
}

void TestIsing2D::check_block_spin_transformation(int n_samples_eq,
                                                  int b,
                                                  int N,
                                                  vec2D K){
    
    Ising2D* pIsing = new Ising2D(N, K);
    pIsing->equilibrate(n_samples_eq);
    pIsing->display_spins();
    
    int n_transformations = floor(log(N)/log(b))-1;
    for(int n = 0; n < n_transformations; ++n){
        pIsing = pIsing->block_spin_transformation(b);
        pIsing->display_spins();
    }
}

#include "test_ising.hpp"

void TestIsing2D::equilibrate_sweepK1(int n_samples,
                                      int N,
                                      double k0,
                                      double kf,
                                      double dk){
    
    double eps = 1E-5;
    for(double k = k0; k <= kf+eps; k += dk){
        
        vec2D K;
        K(0) = k;
        K(1) = 0.0;
        
        Ising2D* pIsing = new Ising2D(N, K);
        pIsing->equilibrate(n_samples);
    }
}

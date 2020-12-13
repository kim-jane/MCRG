#pragma once
#include "ising.hpp"

class MonteCarloRenormalizationGroup{
    
public:
    
    MonteCarloRenormalizationGroup(int b);
    ~MonteCarloRenormalizationGroup(){}
    
    int b_;
    int verbose_;
    int rank_;
    int n_processes_;
    
    double calc_critical_exponent(int n_samples,
                                  int N,
                                  vec2D Kc);

    vec2D locate_critical_point(int n_iterations,
                                int n_samples,
                                int L,
                                vec2D K0);
    
    vec2D approx_critical_point(int n_samples,
                                int L,
                                vec2D K);
    
    int split_samples(int n_samples);
};

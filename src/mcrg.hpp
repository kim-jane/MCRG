#pragma once
#include "ising.hpp"

class MonteCarloRenormalizationGroup{
    
public:
    
    MonteCarloRenormalizationGroup(int b);
    ~MonteCarloRenormalizationGroup(){}
    

    int b_;
    bool verbose_;
    
    void run(int n_samples,
             int N,
             vec2D K);

    vec2D approx_critical_point(int n_samples,
                                int N,
                                vec2D K);
    //void calc_critical_exponent();
    //imat block_spin_transformation(imat spins);
};

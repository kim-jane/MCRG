#pragma once
#include "ising.hpp"

class MonteCarloRenormalizationGroup{
    
public:
    
    MonteCarloRenormalizationGroup(int b,
                                   double T);
    ~MonteCarloRenormalizationGroup(){}
    
    

    
    int b_;
    double T_;

    
    void locate_critical_point(int n_samples,
                               int N0,
                               vec2D K0);
    //void calc_critical_exponent();
    //imat block_spin_transformation(imat spins);
};

#pragma once
#include "ising.hpp"

class TestIsing2D{
    
public:
    
    TestIsing2D(){}
    ~TestIsing2D(){}
    
    void equilibrate_sweepT(int n_samples,
                            int N,
                            double T0,
                            double Tf,
                            double dT);
    
    void check_block_spin_transformation(int n_samples_eq,
                                         int b,
                                         int N,
                                         vec2D K);
};

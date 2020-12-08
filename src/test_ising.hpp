#pragma once
#include "ising.hpp"

class TestIsing2D{
    
public:
    
    TestIsing2D(){}
    ~TestIsing2D(){}
    
    void equilibrate_sweepK1(int n_samples,
                             int N,
                             double k0,
                             double kf,
                             double dk);
};

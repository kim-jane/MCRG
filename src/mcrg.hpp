#pragma once
#include "lattice.hpp"
#include "ising.hpp"
#include <memory>

class MonteCarloRenormalizationGroup{
    
public:
    
    MonteCarloRenormalizationGroup(int b);
    ~MonteCarloRenormalizationGroup(){}
    
    int b_;
    int rank_;
    int n_processes_;
    int iter_;
    FILE* fptr_;
    
    void calc_critical_exponent(int n_samples_eq,
                                int n_samples,
                                int N,
                                vec2D Kc);

    vec2D locate_critical_point(int n_iterations,
                                int n_samples_eq,
                                int n_samples,
                                int L,
                                vec2D K0);
    
    vec2D approx_critical_point(int n_samples_eq,
                                int n_samples,
                                int L,
                                vec2D K);
    
    std::shared_ptr<Lattice> block_spin_transformation(std::shared_ptr<Lattice> pLattice);
    
    int split_samples(int n_samples);
};

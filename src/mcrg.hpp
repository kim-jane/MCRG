#pragma once
#include "lattice.hpp"
#include "ising.hpp"
#include <memory>

class MonteCarloRenormalizationGroup{
    
public:
    
    MonteCarloRenormalizationGroup(int b);
    ~MonteCarloRenormalizationGroup(){}
    
    int b_;

    
    void calc_critical_exponent(int n_samples_eq,
                                int n_samples,
                                int N,
                                double K);

    void locate_critical_point(int n_iterations,
                                int n_samples_eq,
                                int n_samples,
                                int L,
                                double K0);
    
private:
    
    int rank_;
    int n_processes_;
    int iter_;
    FILE* fptr_;
    
    double approx_critical_point(int n_samples_eq,
                                 int n_samples,
                                 int L,
                                 double K);
    
    std::shared_ptr<Lattice> block_spin_transformation(std::shared_ptr<Lattice> pLattice);

    int split_samples(int n_samples);
};

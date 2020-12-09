#pragma once
#include <mpi.h>
#include "definitions.hpp"

class Ising2D{
public:
    
    unsigned n_spins_;
    unsigned N_;
    vec K_;
    imat spins_;
    int rank_;
    int n_processes_;
    int a_;
    
    std::uniform_int_distribution<int> rand_spin_index_;
    
    Ising2D(unsigned N, vec2D K);
    ~Ising2D(){}
    
    void sample_spins();
    void set_spins(imat new_spins);
    void equilibrate(int n_samples, bool write);
    double calc_magnetization();
    double calc_energy();
    vec2D calc_spin_interactions();
    Ising2D* block_spin_transformation(int b);
    void display_spins();
    
private:
    
    void initialize_spins();
    double calc_probability_flip(int i, int j);
    imat nearest_neighbors(int i, int j);
    imat next_nearest_neighbors(int i, int j);
};

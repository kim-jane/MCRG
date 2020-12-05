#pragma once
#include "definitions.hpp"

class Ising2D{
public:
    
    unsigned n_spins_;
    unsigned N_;
    vec K_;
    imat spins_;
    
    std::uniform_int_distribution<int> rand_spin_index_;
    
    Ising2D(unsigned N, vec2D K);
    ~Ising2D(){}
    
    void sample_spins();
    void set_spins(imat new_spins);
    void initialize_spins(FILE* fptr);
    void equilibrate(FILE* fptr);
    double calc_magnetization();
    double calc_energy();
    double calc_probability_flip(int i, int j);
    vec2D calc_correlation();
    imat nearest_neighbors(int i, int j);
    imat next_nearest_neighbors(int i, int j);
    Ising2D* block_spin_transformation(int b);
    void display_spins();
};

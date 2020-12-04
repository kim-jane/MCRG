#pragma once
#include "definitions.hpp"

class Ising2D{
public:
    
    unsigned n_spins_;
    unsigned N_;
    double T_;
    double beta_;
    vec K_;
    imat spins_;
    
    std::uniform_int_distribution<int> rand_spin_index_;
    std::string name_;
    FILE* fptr_;
    
    Ising2D(unsigned N, double T, vec2D K);
    ~Ising2D(){}
    
    void flip_spins();
    void set_spins(imat new_spins);
    void initialize_spins();
    void equilibrate(int n_cycles);
    int calc_magnetization();
    double calc_energy();
    double calc_probability_flip(int i, int j);
    vec2D calc_correlation();
    imat nearest_neighbors(int i, int j);
    imat next_nearest_neighbors(int i, int j);
    Ising2D* block_spin_transformation(int b);
};

#pragma once
#include <mpi.h>
#include "definitions.hpp"

class Ising2D{
public:
    
    int N_;
    int n_spins_;
    int a_;
    vec K_;
    imat spins_;

    std::uniform_int_distribution<int> rand_spin_index_;
    
    Ising2D(int N, vec2D K);
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
    
    int rank_;
    int n_processes_;
    
    void initialize_spins();
    double calc_probability_flip(int i, int j);
    imat nearest_neighbors(int i, int j);
    imat next_nearest_neighbors(int i, int j);
};

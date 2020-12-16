#pragma once
#include <mpi.h>
#include "definitions.hpp"

class Lattice{
    
public:
    
    Lattice(int N);
    Lattice(int a, imat spins);
    ~Lattice(){}
    
    int N_;
    int n_spins_;
    int a_;
    imat spins_;
    
    void write_spins(FILE* fptr);
    void display_spins();
    int choose_random_spin();
    double calc_nearest_neighbor_interaction();
    vec2D calc_interactions();
    imat nearest_neighbors(int i, int j);
    imat next_nearest_neighbors(int i, int j);
    
private:
    
    int rank_;
    int n_processes_;
    double Snn_;
    vec2D S_;
    imat nn_, nnn_;

    std::uniform_int_distribution<int> rand_spin_index_;
    
    void initialize_random_spins();
};

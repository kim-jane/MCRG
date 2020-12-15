#pragma once
#include <mpi.h>
#include "definitions.hpp"
#include "lattice.hpp"

class IsingModel{
public:
    
    IsingModel(vec2D K);
    ~IsingModel(){}
    
    int rank_;
    int n_processes_;
    vec K_;
    
    std::vector<int> cluster_;

    void equilibrate(Lattice* pLattice, int n_samples_eq, bool write);
    double calc_magnetization(Lattice* pLattice);
    double calc_energy(Lattice* pLattice);
    void sample_new_configuration(Lattice* pLattice);
    bool grow_cluster(Lattice* pLattice);
    bool in_cluster(int k);
    

};

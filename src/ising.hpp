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

    void equilibrate(std::shared_ptr<Lattice>, int n_samples_eq, bool write);
    double calc_magnetization(std::shared_ptr<Lattice>);
    double calc_energy(std::shared_ptr<Lattice>);
    void sample_new_configuration(std::shared_ptr<Lattice>);
    bool grow_cluster(std::shared_ptr<Lattice>);
    bool in_cluster(int k);
    

};

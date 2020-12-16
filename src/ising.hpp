#pragma once
#include <mpi.h>
#include <memory>
#include "definitions.hpp"
#include "lattice.hpp"

// 2D Ising Model with nearest neighbor interaction only.
class IsingModel{
    
public:
    
    IsingModel(double K);
    ~IsingModel(){}
    
    void equilibrate(std::shared_ptr<Lattice>,
                     int n_samples_eq,
                     bool write);
    void sample_new_configuration(std::shared_ptr<Lattice> pLattice);

    
private:
    
    int rank_;
    int n_processes_;
    int n_cluster_;
    double K_, P_;
    double E_, E2_, E_avg_, E2_avg_;
    double M_, M2_, M_avg_, M2_avg_;
    double E_sigma_, M_sigma_;
    double C_, Chi_;
    bool cluster_grew_;
    std::vector<int> cluster_;
    imat nn_;
    FILE* fptr_;
    
    void grow_cluster(std::shared_ptr<Lattice> pLattice);
    bool in_cluster(int k);
    double calc_magnetization(std::shared_ptr<Lattice> pLattice);
    double calc_energy(std::shared_ptr<Lattice> pLattice);
};

#pragma once
#include "lattice.hpp"
#include "ising.hpp"
#include <memory>


class RenormalizationGroupNeuralNetwork{
public:
    
    RenormalizationGroupNeuralNetwork(int b);
    ~RenormalizationGroupNeuralNetwork(){}
    
    int n_processes_;
    int rank_;
    int b_;
    int t_;
    double eta_;
    double beta1_;
    double beta2_;
    double epsilon_;
    mat m_;
    mat v_;
    mat W_;
    FILE* fptr_;

    
    // random initial weights and optimizer
    void initialize();
    
    // train RGNN over a range of couplings
    void train_energy(int L,
                      int n_cycles,
                      int n_samples,
                      int n_samples_eq,
                      double T,
                      double h,
                      double eta,
                      double lambda);
    
    // forward-pass
    double scalar_output(const imat& input_spins);
    
    void apply_filter(mat& input);
    
    // simple finite-difference derivative
    mat calc_gradient_scalar_output(double h,
                                  const imat& input_spins);
    
    // simple gradient descent
    void update_weights(double eta,
                        const mat& gradient);

};


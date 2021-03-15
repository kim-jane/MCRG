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
    mat W_;
    
    // train RGNN over a range of couplings
    void train(int N,
               int n_cycles,
               int n_samples,
               int n_samples_eq,
               double T);
    
    // forward-pass
    double predict_temperature(const imat& input_spins);
    
    void apply_filter(mat& input);
    
    // simple finite-difference derivative
    mat calc_temperature_gradient(double h,
                                  const imat& input_spins);
    
    // simple gradient descent
    void update_weights(double eta,
                        const mat& gradient);

};


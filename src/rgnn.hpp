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
    double epsilon_ = 1.0E-8;
    vec m_;
    vec v_;
    mat W_;
    FILE* fptr_;

    
    // random initial weights and optimizer
    void initialize();
    
    // train RGNN over a range of couplings
    void train(int N,
               int n_cycles,
               int n_samples,
               int n_samples_eq,
               double T,
               double h,
               double eta);
    
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


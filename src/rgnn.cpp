#include "rgnn.hpp"

RenormalizationGroupNeuralNetwork::RenormalizationGroupNeuralNetwork(int b){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    b_ = b;
    W_.resize(b,b);
    for(int i = 0; i < b; ++i){
        for(int j = 0; j < b; ++j){
            W_(i,j) = rand_weight();
        }
    }
    
    MPI_Bcast(W_.data(), b_*b_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}



// train RGNN over a range of couplings
void RenormalizationGroupNeuralNetwork::train(int N,
                                              int n_cycles,
                                              int n_samples,
                                              int n_samples_eq,
                                              double T,
                                              double h,
                                              double eta){
    
    int n_samples_loc = split_samples(rank_, n_processes_, n_samples);
    double mse, mse_loc;
    double T_pred, T_pred_avg, T_pred_avg_loc;
    double T_pred_sigma, T_pred2_avg, T_pred2_avg_loc;
    mat grad(b_,b_);
    mat grad_loc(b_,b_);
    
    // initialize Ising model at chosen coupling
    std::unique_ptr<IsingModel> pIsing(new IsingModel(-1.0/T));
    
    // equilibrate lattice
    std::shared_ptr<Lattice> pLattice(new Lattice(N));
    pIsing->equilibrate(pLattice, n_samples_eq, false);
        
    for(int cycles = 0; cycles < n_cycles; ++cycles){
        
        // calculate mse and gradient
        T_pred_avg_loc = 0.0;
        T_pred2_avg_loc = 0.0;
        mse_loc = 0.0;
        grad_loc.setZero();
        for(int samples = 0; samples < n_samples_loc; ++samples){
            
            // get new sample
            pIsing->sample_new_configuration(pLattice);
            
            // pass thru RGNN
            T_pred = predict_temperature(pLattice->spins_);
            T_pred_avg_loc += T_pred;
            T_pred2_avg_loc += T_pred*T_pred;
            
            // calculate
            mse_loc += (T_pred-T)*(T_pred-T);
            grad_loc += 2*(T_pred-T)*calc_temperature_gradient(h, pLattice->spins_);
            
        }
        MPI_Allreduce(&T_pred_avg_loc, &T_pred_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&T_pred2_avg_loc, &T_pred2_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mse_loc, &mse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(grad_loc.data(), grad.data(), b_*b_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // average
        T_pred_avg /= n_samples;
        T_pred2_avg /= n_samples;
        T_pred_sigma = sqrt(T_pred2_avg-T_pred_avg*T_pred_avg);
        mse /= n_samples;
        grad /= n_samples;
        
        if(rank_ == 0){
            printf("%10i %10.5lf %10.5lf %10.5lf %10.5lf\n", cycles, T_pred_avg, T_pred_sigma, mse, grad.norm());
        }
        
        
        // gradient descent
        update_weights(eta, grad);
    }
    
    if(rank_ == 0) std::cout << W_ << std::endl;
}

// forward-pass
double RenormalizationGroupNeuralNetwork::predict_temperature(const imat& input_spins){
    
    mat output = input_spins.cast<double>();
    while(output.rows() > 1){
        apply_filter(output);
    }
    return output(0,0);
}

void RenormalizationGroupNeuralNetwork::apply_filter(mat& input){
    
    int N = input.rows()/b_;
    mat output(N,N);
    mat conv(b_,b_);
    output.setZero();
    
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            
            conv = W_*input.block(b_*i, b_*j, b_, b_);
            
            for(int k = 0; k < b_; ++k){
                for(int l = 0; l < b_; ++l){
                    
                    // relu act func
                    if(conv(k,l) > 0){
                        output(i,j) += conv(k,l);
                    }
                }
            }
        }
    }
    
    input = output;
}

// simple finite-difference derivative
mat RenormalizationGroupNeuralNetwork::calc_temperature_gradient(double h,
                                                                 const imat& input_spins){
    
    double W, T1, T2;
    
    mat grad(b_, b_);
    grad.setZero();
    
    for(int i = 0; i < b_; ++i){
        for(int j = 0; j < b_; ++j){
            
            // store original weight
            W = W_(i,j);
            
            // finite differences O(h^2)
            W_(i,j) = W+h;
            T1 = predict_temperature(input_spins);
            
            W_(i,j) = W-h;
            T2 = predict_temperature(input_spins);
            
            grad(i,j) = (T1-T2)/(2*h);
            
            // return weight to original value
            W_(i,j) = W;
        }
    }
            
    return grad;
}

// simple gradient descent
void RenormalizationGroupNeuralNetwork::update_weights(double eta,
                                                       const mat& gradient){
    W_ -= eta*gradient;
}


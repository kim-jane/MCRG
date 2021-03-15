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
}



// train RGNN over a range of couplings
void RenormalizationGroupNeuralNetwork::train(int N,
                                              int n_cycles,
                                              int n_samples,
                                              int n_samples_eq,
                                              double K){
    
    int n_samples_loc = split_samples(rank_, n_processes_, n_samples);
    double h = 0.001;
    double eta = 0.001;
    double K_pred, mse, mse_loc;
    mat grad(b_,b_);
    mat grad_loc(b_,b_);
    
    double K_pred_avg, K_pred_avg_loc;
    
    // initialize Ising model at chosen coupling
    std::unique_ptr<IsingModel> pIsing(new IsingModel(K));
    
    // equilibrate lattice
    std::shared_ptr<Lattice> pLattice(new Lattice(N));
    pIsing->equilibrate(pLattice, n_samples_eq, false);
    
        
    for(int cycles = 0; cycles < n_cycles; ++cycles){
        
        // pick a random coupling in specified range
        //K = Ki+(Kf-Ki)*rand_unif();

        
        K_pred_avg_loc = 0.0;
        
        // calculate mse and gradient
        mse_loc = 0.0;
        grad_loc.setZero();
        for(int samples = 0; samples < n_samples_loc; ++samples){
            
            // get new sample
            pIsing->sample_new_configuration(pLattice);
            
            // pass thru RGNN
            K_pred = predict_coupling(pLattice->spins_);
            K_pred_avg_loc += K_pred;
            
            // calculate
            mse_loc += (K_pred-K)*(K_pred-K);
            grad_loc += 2*(K_pred-K)*calc_coupling_gradient(h, pLattice->spins_);
            
        }
        MPI_Allreduce(&K_pred_avg_loc, &K_pred_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mse_loc, &mse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(grad_loc.data(), grad.data(), b_*b_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // average
        K_pred_avg /= n_samples;
        mse /= n_samples;
        grad /= n_samples;
        
        if(rank_ == 0){
            std::cout << cycles << "\t" << K_pred_avg << "\t" << mse << "\t" << grad.norm() << std::endl;
        }
        
        
        // gradient descent
        update_weights(eta, grad);
    }
}

// forward-pass
double RenormalizationGroupNeuralNetwork::predict_coupling(const imat& input_spins){
    
    mat output = input_spins.cast<double>();
    while(output.rows() > 1){
        apply_filter(output);
    }
    return -output(0,0);
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
mat RenormalizationGroupNeuralNetwork::calc_coupling_gradient(double h,
                                                              const imat& input_spins){
    
    double W, K1, K2;
    
    mat grad(b_, b_);
    grad.setZero();
    
    for(int i = 0; i < b_; ++i){
        for(int j = 0; j < b_; ++j){
            
            // store original weight
            W = W_(i,j);
            
            // finite differences
            W_(i,j) = W+h;
            K1 = predict_coupling(input_spins);
            
            W_(i,j) = W-h;
            K2 = predict_coupling(input_spins);
            
            grad(i,j) = (K1-K2)/(2*h);
            
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


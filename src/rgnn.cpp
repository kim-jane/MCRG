#include "rgnn.hpp"

// train with energy instead of temperature
// print variance not not stddev
// add L1 (lasso) regularization to drive weights to zero
// when printing final weights, also print out example values of intermediate conv layers so i can make imshow plot

RenormalizationGroupNeuralNetwork::RenormalizationGroupNeuralNetwork(int b){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    b_ = b;
    beta1_ = 0.9;
    beta2_ = 0.999;
    epsilon_ = 1E-8;
    initialize();
}


// initialize random weights and optimizer
void RenormalizationGroupNeuralNetwork::initialize(){
    
    W_.resize(b_, b_);
    for(int i = 0; i < b_; ++i){
        for(int j = 0; j < b_; ++j){
            W_(i,j) = rand_weight();
        }
    }
    
    MPI_Bcast(W_.data(), b_*b_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    t_ = 1;
    m_.resize(b_, b_);
    v_.resize(b_, b_);
    m_.setZero();
    v_.setZero();
}



// train RGNN to predict energy
void RenormalizationGroupNeuralNetwork::train_energy(int L,
                                                     int n_cycles,
                                                     int n_samples,
                                                     int n_samples_eq,
                                                     double K,
                                                     double h,
                                                     double eta){
    
    // open file
    if(rank_ == 0){
        std::string filename = "trainE_b"+std::to_string(b_)
                               +"_N"+std::to_string(L)
                               +"_K"+get_rounded_str(K, 7)+".txt";
        fptr_ = fopen(filename.c_str(), "w");
        fprintf(fptr_, "# Cycles, Average Predicted Energy, Stddev Predicted Energy, Mean Squared Error, || MSE Gradient ||\n");
    }
    
    int n_samples_loc = split_samples(rank_, n_processes_, n_samples);
    double mse, mse_loc;
    double E_pred, E_pred_avg, E_pred_avg_loc;
    double E_pred_sigma, E_pred2_avg, E_pred2_avg_loc;
    mat grad(b_,b_);
    mat grad_loc(b_,b_);
    
    double E;
    
    // initialize Ising model at chosen coupling
    std::unique_ptr<IsingModel> pIsing(new IsingModel(K));
    
    // equilibrate large lattice
    std::shared_ptr<Lattice> pLattice(new Lattice(L));
    pIsing->equilibrate(pLattice, n_samples_eq, false);
    
    // equilibrate small lattice
    /*
    int S = L/b_;
    std::shared_ptr<Lattice> pLattice(new Lattice(S));
    pIsing->equilibrate(pLattice, n_samples_eq, false);
     */
        
    for(int cycles = 0; cycles <= n_cycles; ++cycles){
        
        // calculate mse and gradient
        E_pred_avg_loc = 0.0;
        E_pred2_avg_loc = 0.0;
        mse_loc = 0.0;
        grad_loc.setZero();
        for(int samples = 0; samples < n_samples_loc; ++samples){
            
            // get new sample
            pIsing->sample_new_configuration(pLattice);
            
            // calculate nn energy per spin
            E = pIsing->calc_energy(pLattice);
            
            // pass sample thru RGNN
            E_pred = scalar_output(pLattice->spins_);
            E_pred_avg_loc += E_pred;
            E_pred2_avg_loc += E_pred*E_pred;
            
            // calculate cost and gradient
            mse_loc += (E_pred-E)*(E_pred-E);
            grad_loc += 2*(E_pred-E)*calc_gradient_scalar_output(h, pLattice->spins_);
        
        }
        MPI_Allreduce(&E_pred_avg_loc, &
                      E_pred_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&E_pred2_avg_loc, &E_pred2_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mse_loc, &mse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(grad_loc.data(), grad.data(), b_*b_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // average
        E_pred_avg /= n_samples;
        E_pred2_avg /= n_samples;
        E_pred_sigma = sqrt(E_pred2_avg-E_pred_avg*E_pred_avg);
        mse /= n_samples;
        grad /= n_samples;
        
        if(rank_ == 0){
            fprintf(fptr_, "%10i%15.7e%15.7e%15.7e%15.7e\n", cycles, E_pred_avg, E_pred_sigma, mse, grad.norm());
        }
        
        
        // gradient descent
        update_weights(eta, grad);
    }
    
    if(rank_ == 0){
        
        fprintf(fptr_, "\nFinal Weights: ");
        for(int i = 0; i < b_; ++i){
            for(int j = 0; j < b_; ++j){
                fprintf(fptr_, "%20.10lf", W_(i,j));
            }
        }
        fclose(fptr_);
    }
}

// recursively apply filter until only one scalar remains
double RenormalizationGroupNeuralNetwork::scalar_output(const imat& input_spins){
    
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
mat RenormalizationGroupNeuralNetwork::calc_gradient_scalar_output(double h,
                                                                   const imat& input_spins){
    
    double W, out1, out2;
    
    mat grad(b_, b_);
    grad.setZero();
    
    for(int i = 0; i < b_; ++i){
        for(int j = 0; j < b_; ++j){
            
            // store original weight
            W = W_(i,j);
            
            // finite differences O(h^2)
            W_(i,j) = W+h;
            out1 = scalar_output(input_spins);
            
            W_(i,j) = W-h;
            out2 = scalar_output(input_spins);
            
            grad(i,j) = (out1-out2)/(2*h);
            
            // return weight to original value
            W_(i,j) = W;
        }
    }
            
    return grad;
}

// simple gradient descent
void RenormalizationGroupNeuralNetwork::update_weights(double eta,
                                                       const mat& gradient){
    
    eta_ = eta*sqrt(1.0-pow(beta2_, t_))/(1.0-pow(beta1_, t_));
    
    for(int i = 0; i < b_; ++i){
        for(int j = 0; j < b_; ++j){
            m_(i,j) = beta1_*m_(i,j)+(1.0-beta1_)*gradient(i,j);
            v_(i,j) = beta2_*v_(i,j)+(1.0-beta2_)*gradient(i,j)*gradient(i,j);
            W_(i,j) -= eta_*m_(i,j)/(sqrt(v_(i,j))+epsilon_);
        }
    }
    
    t_ += 1;
}


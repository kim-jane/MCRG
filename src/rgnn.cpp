#include "rgnn.hpp"



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

// train so that scalar output of network is scale invariant

void RenormalizationGroupNeuralNetwork::train_scalar_output(int L,
                                                            int n_cycles,
                                                            int n_samples,
                                                            int n_samples_eq,
                                                            double K,
                                                            double h,
                                                            double eta){
    
    // open file
    if(rank_ == 0){
        std::string filename = "train_scalar_b"+std::to_string(b_)
                               +"_L"+std::to_string(L)
                               +"_K"+get_rounded_str(K, 7)+".txt";
        fptr_ = fopen(filename.c_str(), "w");
        
        fprintf(fptr_, "# Initial Weights: ");
        for(int i = 0; i < b_; ++i){
            for(int j = 0; j < b_; ++j){
                fprintf(fptr_, "%20.10lf", W_(i,j));
            }
        }
        
        
        fprintf(fptr_, "\n# Cycles, Avg Output L, Var Output L, Avg Output S, Var Output S, MSE, || MSE Gradient ||\n");
        

    }

    int n_samples_loc = split_samples(rank_, n_processes_, n_samples);
    double uL, uL_avg, uL_avg_loc, uL2_avg, uL2_avg_loc, uL_var;
    double uS, uS_avg, uS_avg_loc, uS2_avg, uS2_avg_loc, uS_var;
    double mse;
    mat gradL(b_,b_);
    mat gradL_loc(b_,b_);
    mat gradS(b_,b_);
    mat gradS_loc(b_,b_);
    mat grad(b_,b_);
    
    // initialize Ising model at chosen coupling
    std::unique_ptr<IsingModel> pIsing(new IsingModel(K));
    
    // equilibrate large lattices
    std::shared_ptr<Lattice> pLatticeL(new Lattice(L));
    pIsing->equilibrate(pLatticeL, n_samples_eq, false);
    
    // equilibrate small lattice
    int S = L/b_;
    std::shared_ptr<Lattice> pLatticeS(new Lattice(S));
    pIsing->equilibrate(pLatticeS, n_samples_eq, false);
    
    // minimize
    for(int cycles = 0; cycles <= n_cycles; ++cycles){
        
        uL_avg_loc = 0.0;
        uL2_avg_loc = 0.0;
        uS_avg_loc = 0.0;
        uS2_avg_loc = 0.0;
        gradL_loc.setZero();
        gradS_loc.setZero();
        
        for(int samples = 0; samples < n_samples_loc; ++samples){
            
            // get new sample
            pIsing->sample_new_configuration(pLatticeL);
            pIsing->sample_new_configuration(pLatticeS);
            
            // pass sample thru RGNN
            uL = scalar_output(pLatticeL->spins_);
            uS = scalar_output(pLatticeS->spins_);
            
            // add up values for averages
            uL_avg_loc += uL;
            uS_avg_loc += uS;
            uL2_avg_loc += uL*uL;
            uS2_avg_loc += uS*uS;
            gradL_loc += calc_gradient_scalar_output(h, pLatticeL->spins_);
            gradS_loc += calc_gradient_scalar_output(h, pLatticeS->spins_);
        }
        
        MPI_Allreduce(&uL_avg_loc, &uL_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&uS_avg_loc, &uS_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&uL2_avg_loc, &uL2_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&uS2_avg_loc, &uS2_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(gradL_loc.data(), gradL.data(), b_*b_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(gradS_loc.data(), gradS.data(), b_*b_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // finalize averages
        uL_avg /= n_samples;
        uS_avg /= n_samples;
        uL2_avg /= n_samples;
        uS2_avg /= n_samples;
        gradL /= n_samples;
        gradS /= n_samples;
        
        // variances
        uL_var = uL2_avg-uL_avg*uL_avg;
        uS_var = uS2_avg-uS_avg*uS_avg;
        
        // cost and gradient
        mse = (uL_avg-uS_avg)*(uL_avg-uS_avg);
        grad = 2*(uL_avg-uS_avg)*(gradL-gradS);
        
        
        if(rank_ == 0){
            if(cycles%10 == 0){
                printf("%10i%15.7e%15.7e%15.7e%15.7e\n", cycles, uL_avg, uS_avg, mse, grad.norm());
            }
            
            fprintf(fptr_, "%10i%15.7e%15.7e%15.7e%15.7e%15.7e%15.7e\n", cycles, uL_avg, uL_var, uS_avg, uS_var, mse, grad.norm());
        }
        
        
        // gradient descent
        update_weights(eta, grad);
    }
    
    if(rank_ == 0){
        
        fprintf(fptr_, "\n#Final Weights: ");
        for(int i = 0; i < b_; ++i){
            for(int j = 0; j < b_; ++j){
                fprintf(fptr_, "%20.10lf", W_(i,j));
            }
        }
        fclose(fptr_);
    }
}


void RenormalizationGroupNeuralNetwork::test_temperature(int N,
                                                         int n_samples,
                                                         int n_samples_eq,
                                                         double K,
                                                         double DeltaK){
    
    // open file
    if(rank_ == 0){
        std::string filename = "testT_b"+std::to_string(b_)
                               +"_N"+std::to_string(N)
                               +"_K"+get_rounded_str(K, 7)+".txt";
        fptr_ = fopen(filename.c_str(), "w");
        fprintf(fptr_, "# Temperature, Average Predicted Temperature, Variance Predicted Temperature, Mean Squared Error\n");
    }
    
    int n_samples_loc = split_samples(rank_, n_processes_, n_samples);
    double mse, mse_loc;
    double T_pred, T_pred_avg, T_pred_avg_loc;
    double T_pred_var, T_pred2_avg, T_pred2_avg_loc;
    double dK = 2*DeltaK/15;
    
    for(double k = K-DeltaK; k <= K+DeltaK; k += dK){
        
        // corresponding temperature
        double T = -1.0/k;
        
        // initialize Ising model at chosen coupling
        std::unique_ptr<IsingModel> pIsing(new IsingModel(k));
        
        // equilibrate large lattice
        std::shared_ptr<Lattice> pLattice(new Lattice(N));
        pIsing->equilibrate(pLattice, n_samples_eq, false);
        
        // calculate mse and gradient
        T_pred_avg_loc = 0.0;
        T_pred2_avg_loc = 0.0;
        mse_loc = 0.0;
        for(int samples = 0; samples < n_samples_loc; ++samples){
            
            // get new sample
            pIsing->sample_new_configuration(pLattice);
            
            // pass sample thru RGNN
            T_pred = scalar_output(pLattice->spins_);
            T_pred_avg_loc += T_pred;
            T_pred2_avg_loc += T_pred*T_pred;
            
            // calculate cost and gradient
            mse_loc += (T_pred-T)*(T_pred-T);
        
        }
        MPI_Allreduce(&T_pred_avg_loc, &
                      T_pred_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&T_pred2_avg_loc, &T_pred2_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mse_loc, &mse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // average
        T_pred_avg /= n_samples;
        T_pred2_avg /= n_samples;
        T_pred_var = T_pred2_avg-T_pred_avg*T_pred_avg;
        mse /= n_samples;
        
        
        if(rank_ == 0){
            std::cout << T << "\t" << T_pred_avg << "\t" << mse << std::endl;
            fprintf(fptr_, "%15.7e%15.7e%15.7e%15.7e\n", T, T_pred_avg, T_pred_var, mse);
        }
    }
    
    if(rank_ == 0) fclose(fptr_);
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
            output(i,j) = conv.maxCoeff();
            //output(i,j) = exp(-conv.sum());

            /*
            for(int k = 0; k < b_; ++k){
                for(int l = 0; l < b_; ++l){

                    // relu act func and sum
                    if(conv(k,l) > 0){
                        output(i,j) += conv(k,l);
                    }
                }
            }
             */
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
    
    MPI_Bcast(W_.data(), b_*b_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    t_ += 1;
}


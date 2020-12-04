#include "mcrg.hpp"

MonteCarloRenormalizationGroup::MonteCarloRenormalizationGroup(int b,
                                                               double T){
    
    b_ = b;
    T_ = T;
}


void MonteCarloRenormalizationGroup::locate_critical_point(int n_samples,
                                                           int N,
                                                           vec2D K){
    
    if(N%b_ != 0){
        printf("Lattice size is not divisible by the scaling factor.\n");
        exit(1);
    }

    
    // two lattices with sizes differing by b
    Ising2D* pIsing1 = new Ising2D(N, T_, K);
    Ising2D* pIsing2 = new Ising2D(N/b_, T_, K);
    
    // equilibrate both lattices
    int n_cycles_eq = 1E5;
    pIsing1->equilibrate(n_cycles_eq);
    pIsing2->equilibrate(n_cycles_eq);
    
    // transform larger lattice
    Ising2D* pIsing1b = pIsing1->block_spin_transformation(b_);
    
    // calculate correlation functions
    vec2D S1, S2, S1b;
    vec2D S1_avg, S2_avg, S1b_avg;
    mat2D S1b_S1_avg, S2_S2_avg;
    S1_avg.setZero();
    S2_avg.setZero();
    S1b_avg.setZero();
    S1b_S1_avg.setZero();
    S2_S2_avg.setZero();
    
    for(int samples = 0; samples < n_samples; ++samples){
        
        // get new samples
        pIsing1->flip_spins();
        pIsing2->flip_spins();
        pIsing1b->flip_spins();
        
        // calculate correlation for each lattice
        S1 = pIsing1->calc_correlation();
        S2 = pIsing2->calc_correlation();
        S1b = pIsing1b->calc_correlation();
        
        // add up values for averages
        S1_avg += S1/n_samples;
        S2_avg += S2/n_samples;
        S1b_avg += S1b/n_samples;
        S1b_S1_avg += S1b*S1.transpose()/n_samples;
        S2_S2_avg += S2*S2.transpose()/n_samples;
    }
    
    // divide by 4?
    
    /*
    S1_avg /= n_samples;
    S2_avg /= n_samples;
    S1b_avg /= n_samples;
    S1b_S1_avg /= n_samples;
    S2_S2_avg /= n_samples;
    */
    
    mat2D dS1_dK = S1b_S1_avg-S1b_avg*S1_avg.transpose();
    mat2D dS2_dK = S2_S2_avg-S2_avg*S2_avg.transpose();
    
    vec2D dK = (dS1_dK-dS2_dK).inverse() * (S1b_avg-S2_avg);
    
    printf("%.10lf %.10lf\n", dK(0), dK(1));
}

/*
void MonteCarloRenormalizationGroup::calc_critical_exponent(int n_samples){
    
    std::cout << n_max_ << std::endl;
    std::cout << IsingModels_[0]->spins_ << std::endl;

    for(int n = 1; n < n_max_; ++n){
        
        // RG transformation
        block_spin_transformation(IsingModels_[n-1]->spins_);
        
        // set spins of scaled system
        N_(n) = N_(n-1)/b_;
        J_(n) = J_(n-1);
        K_(n) = K_(n-1);
        T_(n) = T_(n-1);
        IsingModels_.push_back(new Ising2D(N_(n), J_(n), K_(n), T_(n)));
        IsingModels_[n]->set_spins(block_spins_);
        
        // calculate RG transformation matrix
        calc_transformation_matrix(n, n_samples);
        
        std::cout << block_spins_ << std::endl;
    }
}

void MonteCarloRenormalizationGroup::block_spin_transformation(imat spins){
    

}

void MonteCarloRenormalizationGroup::calc_transformation_matrix(int n, int n_samples){
    
    S1_avg_ = 0.0;
    S2_avg_ = 0.0;
    S1_prev_avg_ = 0.0;
    S2_prev_avg_ = 0.0;
    S1_S1_avg_ = 0.0;
    S1_S2_avg_ = 0.0; // same as S2_S1_avg_
    S2_S1_avg_ = 0.0;
    S1_S1_prev_avg_ = 0.0;
    S1_S2_prev_avg_ = 0.0;
    S2_S1_prev_avg_ = 0.0;
    S2_S2_prev_avg_ = 0.0;
    
    
    for(int samples = 0; samples < n_samples; ++samples){
        
        // get new spin configurations
        IsingModels_[n]->flip_spins();
        IsingModels_[n-1]->flip_spins();
        
        // calculate spin interaction terms
        S1_ = 0.0;
        S2_ = 0.0;
        for(int i = 0; i < IsingModels_[n]->N_; ++i){
            for(int j = 0; j < IsingModels_[n]->N_; ++j){
                
                get_nearest_neighbors();
                get_next_nearest_neighbors();
                
                for(int n = 0; n < 4; ++n){
                    E_ += J_*spins_(i_,j_)*spins_(NN_(n,0),NN_(n,1));
                    E_ += K_*spins_(i_,j_)*spins_(NNN_(n,0),NNN_(n,1));
                }
            }
        }
        
        S1_prev_ = 0.0;
        S2_prev_ = 0.0;
        
        
    }
}
*/

/*     if(N0%b == 0){

    n_max_ = floor(log(N0)/log(b_));
    IsingModels_.reserve(n_max_);
    N_.resize(n_max_);
    J_.resize(n_max_);
    K_.resize(n_max_);
    T_.resize(n_max_);
    N_[0] = N0;
    J_[0] = J0;
    K_[0] = K0;
    T_[0] = T0;
    
    // equilibrate original system
    int n_cycles_eq = 1E5;
    IsingModels_.push_back(new Ising2D(N0, J0, K0, T0));
    
    
    printf("Equilibrating initial system...\n");
    IsingModels_[0]->equilibrate(n_cycles_eq);
    printf("Done.\n");
}
else{
    printf("Initial lattice size is not divisible by the scaling factor.\n");
    exit(1);
}
 */

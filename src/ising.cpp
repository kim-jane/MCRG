#include "ising.hpp"

Ising2D::Ising2D(unsigned N,
                 vec2D K){
    
    N_ = N;
    K_ = K;
    
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0,N_-1);
}

void Ising2D::equilibrate(){

    initialize_spins();
    //printf("Initial random spins for N = %i, K = (%.3f, %.3f)\n\n", N_, K_(0), K_(1));
    //display_spins();

    printf("Equilibrating...\n\n");
    printf("%15s %10s %10s \n", "Iteration", "E/spin", "M/spin");
    fflush(stdout);
    
    int i = 0;
    int patience0 = 2000;
    int patience = patience0;
    double E = 1.0;
    double E_min = 1.0;
    double M = 0.0;
    
    while(patience > 0){
        
        sample_spins();
        E = calc_energy();
        M = calc_magnetization();
        
        if(print(i)) printf("%15i %10lf %10lf \n", i, E, M);
        
        if(E < E_min){
            E_min = E;
            patience = patience0;
        }
        else{
            patience--;
        }
            
        i++;
    }
    
    printf("%15i %10lf %10lf \n\n", i, E, M);
    //printf("Equilibrated spins\n\n");
    //display_spins();
}



double Ising2D::calc_energy(){
    
    double E = 0.0;
    imat NN, NNN;
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            NN = nearest_neighbors(i,j);
            NNN = next_nearest_neighbors(i,j);
            
            for(int n = 0; n < 4; ++n){
                E += K_[0]*spins_(i,j)*spins_(NN(n,0),NN(n,1));
                E += K_[1]*spins_(i,j)*spins_(NNN(n,0),NNN(n,1));
            }
        }
    }
    
    return E/(4.0*n_spins_);
}


double Ising2D::calc_magnetization(){
    
    return (double)spins_.sum()/n_spins_;
}

vec2D Ising2D::calc_correlation(){
    
    imat NN, NNN;
    vec2D S;
    S.setZero();
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            NN = nearest_neighbors(i,j);
            NNN = next_nearest_neighbors(i,j);
            
            for(int n = 0; n < 4; ++n){
                S(0) += spins_(i,j)*spins_(NN(n,0),NN(n,1));
                S(1) += spins_(i,j)*spins_(NNN(n,0),NNN(n,1));
            }
        }
    }
    
    return S/4.0;
}

void Ising2D::sample_spins(){
    
    int i, j;
    double P;
    
    for(int n = 0; n < n_spins_; ++n){
        
        i = rand_spin_index_(rng);
        j = rand_spin_index_(rng);
        P = calc_probability_flip(i,j);
        
        if(rand_unif() < P){
            spins_(i,j) *= -1;
        }
    }
}

double Ising2D::calc_probability_flip(int i, int j){
    
    imat NN = nearest_neighbors(i,j);
    imat NNN = next_nearest_neighbors(i,j);
    
    double dE = 0.0;
    for(int n = 0; n < 4; ++n){
        dE += K_[0]*spins_(NN(n,0),NN(n,1));
        dE += K_[1]*spins_(NNN(n,0),NNN(n,1));
    }
    
    dE *= 2.0*spins_(i,j);
    
    return exp(dE);
}


void Ising2D::set_spins(imat new_spins){
    
    spins_ = new_spins;
}


void Ising2D::initialize_spins(){
    
    spins_.resize(N_,N_);
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            spins_(i,j) = rand_spin();
        }
    }
}


imat Ising2D::nearest_neighbors(int i, int j){
    
    imat NN(4,2);
    
    NN(0,0) = (i+1)%N_;
    NN(0,1) = j;
    NN(1,0) = (i-1)%N_;
    NN(1,1) = j;
    NN(2,0) = i;
    NN(2,1) = (j+1)%N_;
    NN(3,0) = i;
    NN(3,1) = (j-1)%N_;
    
    return NN;
}


imat Ising2D::next_nearest_neighbors(int i, int j){
    
    imat NNN(4,2);
    
    NNN(0,0) = (i+1)%N_;
    NNN(0,1) = (j+1)%N_;
    NNN(1,0) = (i-1)%N_;
    NNN(1,1) = (j+1)%N_;
    NNN(2,0) = (i+1)%N_;
    NNN(2,1) = (j-1)%N_;
    NNN(3,0) = (i-1)%N_;
    NNN(3,1) = (j-1)%N_;
    
    return NNN;
}

Ising2D* Ising2D::block_spin_transformation(int b){
    
    int Nb = N_/b;
    int spin_tot;
    imat block_spins(Nb, Nb);
    
    // loop thru blocks
    for(int ib = 0; ib < Nb; ++ib){
        for(int jb = 0; jb < Nb; ++jb){
            
            // loop thru spins in each block
            spin_tot = 0;
            for(int i = ib*b; i < (ib+1)*b; ++i){
                for(int j = jb*b; j < (jb+1)*b; ++j){
                    spin_tot += spins_(i,j);
                }
            }
            
            // majority rule
            if(spin_tot == 0){
                block_spins(ib, jb) = rand_spin();
            }
            else if(spin_tot > 0){
                block_spins(ib, jb) = 1;
            }
            else{
                block_spins(ib, jb) = -1;
            }
        }
    }
    
    Ising2D* pIsing = new Ising2D(Nb, K_);
    pIsing->set_spins(block_spins);
    
    printf("Block spin transformation N = %i --> %i\n\n", N_, Nb);
    pIsing->display_spins();
    
    return pIsing;
}




void Ising2D::display_spins(){
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            if(spins_(i,j) == 1){
                display_spin_up();
            }
            else{
                display_spin_down();
            }
        }
        printf("\n");
    }
    printf("\n");
}

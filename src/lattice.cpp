#include "lattice.hpp"

Lattice::Lattice(int N){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    a_ = 1;
    N_ = N;
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0, n_spins_-1);
    initialize_random_spins();
}

Lattice::Lattice(int a, imat spins){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    a_ = a;
    N_ = spins_.rows();
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0, n_spins_-1);
    spins_ = spins;
}

void Lattice::initialize_random_spins(){
    
    spins_.resize(N_,N_);
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            spins_(i,j) = rand_spin();
        }
    }
}


void Lattice::display_spins(){
    
    if(rank_ == 0){
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
        printf("(Lattice spacing a = %i)\n\n", a_);
    }
}

int Lattice::choose_random_spin(){
    
    return rand_spin_index_(rng);
}

double Lattice::calc_nearest_neighbor_interaction(){
    
    double Snn = 0;
    imat nn;
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            nn = nearest_neighbors(i,j);
            
            for(int k = 0; k < 4; ++k){
                Snn += spins_(i,j)*spins_(nn(k,0),nn(k,1));
            }
        }
    }
    
    return Snn;
}


vec2D Lattice::calc_interactions(){
    
    vec2D S;
    S.setZero();
    imat nn, nnn;
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            nn = nearest_neighbors(i,j);
            nnn = next_nearest_neighbors(i,j);
            
            for(int k = 0; k < 4; ++k){
                S(0) += spins_(i,j)*spins_(nn(k,0),nn(k,1));
                S(1) += spins_(i,j)*spins_(nnn(k,0),nnn(k,1));
            }
        }
    }
    
    return S;
}



imat Lattice::nearest_neighbors(int i, int j){
    
    imat nn(4,2);

    nn(0,0) = ((i+1)%N_+N_)%N_;
    nn(0,1) = j;
    nn(1,0) = ((i-1)%N_+N_)%N_;
    nn(1,1) = j;
    nn(2,0) = i;
    nn(2,1) = ((j+1)%N_+N_)%N_;
    nn(3,0) = i;
    nn(3,1) = ((j-1)%N_+N_)%N_;
    
    return nn;
}


imat Lattice::next_nearest_neighbors(int i, int j){
    
    imat nnn(4,2);
    
    nnn(0,0) = ((i+1)%N_+N_)%N_;
    nnn(0,1) = ((j+1)%N_+N_)%N_;
    nnn(1,0) = ((i-1)%N_+N_)%N_;
    nnn(1,1) = ((j+1)%N_+N_)%N_;
    nnn(2,0) = ((i+1)%N_+N_)%N_;
    nnn(2,1) = ((j-1)%N_+N_)%N_;
    nnn(3,0) = ((i-1)%N_+N_)%N_;
    nnn(3,1) = ((j-1)%N_+N_)%N_;
    
    return nnn;
}


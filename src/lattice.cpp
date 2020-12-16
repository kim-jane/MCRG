#include "lattice.hpp"

Lattice::Lattice(int N){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    a_ = 1;
    N_ = N;
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0, n_spins_-1);
    initialize_random_spins();
    nn_.resize(4,2);
    nnn_.resize(4,2);
}


Lattice::Lattice(int a, imat spins){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    spins_ = spins;
    a_ = a;
    N_ = spins_.rows();
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0, n_spins_-1);
    nn_.resize(4,2);
    nnn_.resize(4,2);
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


void Lattice::write_spins(FILE* fptr){
    
    if(rank_ == 0){
        fprintf(fptr, "\n");
        for(int i = 0; i < N_; ++i){
            fprintf(fptr, "# ");
            for(int j = 0; j < N_; ++j){
                fprintf(fptr, "%i, ", spins_(i,j));
            }
            fprintf(fptr, "\n");
        }
    }
}


int Lattice::choose_random_spin(){
    
    return rand_spin_index_(rng);
}


double Lattice::calc_nearest_neighbor_interaction(){
    
    Snn_ = 0.0;
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            nn_ = nearest_neighbors(i,j);
            
            for(int k = 0; k < 4; ++k){
                Snn_ += spins_(i,j)*spins_(nn_(k,0),nn_(k,1));
            }
        }
    }
    
    return Snn_;
}


vec2D Lattice::calc_interactions(){
    
    S_.setZero();
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            nn_ = nearest_neighbors(i,j);
            nnn_ = next_nearest_neighbors(i,j);
            
            for(int k = 0; k < 4; ++k){
                S_(0) += spins_(i,j)*spins_(nn_(k,0),nn_(k,1));
                S_(1) += spins_(i,j)*spins_(nnn_(k,0),nnn_(k,1));
            }
        }
    }
    
    return S_/4.0;
}



imat Lattice::nearest_neighbors(int i, int j){
    
    nn_(0,0) = ((i+1)%N_+N_)%N_;
    nn_(0,1) = j;
    nn_(1,0) = ((i-1)%N_+N_)%N_;
    nn_(1,1) = j;
    nn_(2,0) = i;
    nn_(2,1) = ((j+1)%N_+N_)%N_;
    nn_(3,0) = i;
    nn_(3,1) = ((j-1)%N_+N_)%N_;
    
    return nn_;
}


imat Lattice::next_nearest_neighbors(int i, int j){
    
    nnn_(0,0) = ((i+1)%N_+N_)%N_;
    nnn_(0,1) = ((j+1)%N_+N_)%N_;
    nnn_(1,0) = ((i-1)%N_+N_)%N_;
    nnn_(1,1) = ((j+1)%N_+N_)%N_;
    nnn_(2,0) = ((i+1)%N_+N_)%N_;
    nnn_(2,1) = ((j-1)%N_+N_)%N_;
    nnn_(3,0) = ((i-1)%N_+N_)%N_;
    nnn_(3,1) = ((j-1)%N_+N_)%N_;
    
    return nnn_;
}


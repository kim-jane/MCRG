#include "ising.hpp"

Ising2D::Ising2D(unsigned N,
                 vec2D K,
                 int verbose){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    N_ = N;
    K_ = K;
    verbose_ = verbose;
    
    a_ = 1;
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0,N_-1);
}

void Ising2D::equilibrate(int n_samples, FILE* fptr){

    initialize_spins();

    if(rank_ == 0){
        
        if(fptr != NULL){
            fprintf(fptr, "# %13s, %10s, %10s, %10s, %10s\n", "Iteration", "E", "E/spin", "M", "M/spin");
        }

        if(verbose_ > 0){
            printf("Equilibrating %i N = %i lattice(s) at K = ", n_processes_, N_);
            print_vec2D(K_);
        }
        if(verbose_ > 1){
            printf("%15s %10s %10s \n", "Iteration", "E/spin", "M/spin");
        }
    }
    
    double E, M, e, m;
    for(int n = 1; n <= n_samples; ++n){
        
        sample_spins();
        E = calc_energy();
        M = calc_magnetization();
        e = E/n_spins_;
        m = M/n_spins_;
        
        if(rank_ == 0){
            
            if(fptr != NULL){
                fprintf(fptr, "%15i, %10.7e, %10.7e, %10.7e, %10.7e,\n", n, E, e, M, m);
            }
            
            if(verbose_ > 1 && print_iter(n)){
                printf("%15i %10lf %10lf \n", n, e, m);
            }
        }
    }

    if(rank_ == 0 && verbose_ > 1){
        printf("\nEquilibrated spins\n");
        display_spins();
    }
}



double Ising2D::calc_energy(){
    
    double E = 0.0;
    imat nn, nnn;
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            nn = nearest_neighbors(i,j);
            nnn = next_nearest_neighbors(i,j);
            
            for(int n = 0; n < 4; ++n){
                E += K_[0]*spins_(i,j)*spins_(nn(n,0),nn(n,1));
                E += K_[1]*spins_(i,j)*spins_(nnn(n,0),nnn(n,1));
            }
        }
    }
    
    return E/4.0;
}


double Ising2D::calc_magnetization(){
    
    return spins_.sum();
}

vec2D Ising2D::calc_spin_interactions(){
    
    imat nn, nnn;
    vec2D S;
    S.setZero();
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            nn = nearest_neighbors(i,j);
            nnn = next_nearest_neighbors(i,j);
            
            for(int n = 0; n < 4; ++n){
                S(0) += spins_(i,j)*spins_(nn(n,0),nn(n,1));
                S(1) += spins_(i,j)*spins_(nnn(n,0),nnn(n,1));
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
    
    imat nn = nearest_neighbors(i,j);
    imat nnn = next_nearest_neighbors(i,j);
    
    double dE = 0.0;
    for(int n = 0; n < 4; ++n){
        dE += K_[0]*spins_(nn(n,0),nn(n,1));
        dE += K_[1]*spins_(nnn(n,0),nnn(n,1));
    }
    
    dE *= 2.0*spins_(i,j);
    
    return exp(dE);
}


void Ising2D::initialize_spins(){
    
    spins_.resize(N_,N_);
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            spins_(i,j) = rand_spin();
        }
    }
    
    if(rank_ == 0 && verbose_ > 1){
        printf("Initial random spins for N = %i\n", N_);
        display_spins();
    }
}


imat Ising2D::nearest_neighbors(int i, int j){
    
    imat nn(4,2);
    
    nn(0,0) = (i+1)%N_;
    nn(0,1) = j;
    nn(1,0) = (i-1)%N_;
    nn(1,1) = j;
    nn(2,0) = i;
    nn(2,1) = (j+1)%N_;
    nn(3,0) = i;
    nn(3,1) = (j-1)%N_;
    
    return nn;
}


imat Ising2D::next_nearest_neighbors(int i, int j){
    
    imat nnn(4,2);
    
    nnn(0,0) = (i+1)%N_;
    nnn(0,1) = (j+1)%N_;
    nnn(1,0) = (i-1)%N_;
    nnn(1,1) = (j+1)%N_;
    nnn(2,0) = (i+1)%N_;
    nnn(2,1) = (j-1)%N_;
    nnn(3,0) = (i-1)%N_;
    nnn(3,1) = (j-1)%N_;
    
    return nnn;
}

Ising2D* Ising2D::block_spin_transformation(int b){
    
    if(N_%b != 0){
        if(rank_ == 0){
            print_error("Lattice size is not divisible by the scaling factor.\n");
            
        }
        exit(1);
    }
    
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
    
    Ising2D* pIsing = new Ising2D(Nb, K_, verbose_);
    pIsing->spins_ = block_spins;
    pIsing->a_ = a_*b;
    
    if(rank_ == 0){
        if(verbose_ > 0){
            printf("Block spin transformation N = %i --> %i\n", N_, Nb);
        }
        if(verbose_ > 1){
            pIsing->display_spins();
        }
    }

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
    printf("Number of spins = %i\n", N_*N_);
    printf("Lattice spacing = %i\n", a_);
}

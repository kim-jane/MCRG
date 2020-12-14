#include "ising.hpp"

Ising2D::Ising2D(int N,
                 vec2D K){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    N_ = N;
    K_ = K;
    a_ = 1;
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0, n_spins_-1);
}

Ising2D::Ising2D(int N,
                 vec2D K){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    N_ = N;
    K_ = K;
    a_ = 1;
    n_spins_ = N_*N_;
    rand_spin_index_ = std::uniform_int_distribution<int>(0, n_spins_-1);
}

void Ising2D::equilibrate(int n_samples_eq){

    FILE* fptr = NULL;
    if(rank_ == 0){
        
        printf("Equilibrating %i N = %i lattice(s) at K = ", n_processes_, N_);
        print_vec2D(K_);
        
        std::string filename = "equilibrate_N_"+std::to_string(N_)
                               +"_K1_"+get_rounded_str(K_(0))
                               +"_K2_"+get_rounded_str(K_(1))
                               +".txt";

        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# Using %i parallel processes\n", n_processes_);
        fprintf(fptr, "# %s, %s, %s, %s, %s, %s, %s\n",
                "Iteration", "Avg E/spin", "Stddev E/spin", "Heat Capacity",
                "Avg |M|/spin", "Stddev |M|/spin", "Susceptibility");
    }

    initialize_spins();
    double E, M, E2, M2;
    double E_avg, M_avg, E2_avg, M2_avg;
    double E_sigma, M_sigma;
    double C, Chi;
    
    for(int n = 0; n <= n_samples_eq; ++n){

        sample_spins();
        
        if(write_iter(n)){
            
            E = calc_energy();
            M = abs(calc_magnetization());
            E2 = E*E;
            M2 = M*M;
            
            MPI_Reduce(&E, &E_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&M, &M_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&E2, &E2_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&M2, &M2_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if(rank_ == 0){
                
                E_avg /= n_processes_;
                M_avg /= n_processes_;
                E2_avg /= n_processes_;
                M2_avg /= n_processes_;
                
                E_sigma = E2_avg-E_avg*E_avg;
                M_sigma = M2_avg-M_avg*M_avg;
                
                C = E_sigma*K_(0)*K_(0);
                Chi = -M_sigma*K_(0);
                
                E_sigma = sqrt(E_sigma);
                M_sigma = sqrt(M_sigma);
                
                fprintf(fptr, "%i, %10.7e, %10.7e, %10.7e, %10.7e, %10.7e, %10.7e\n", n, E_avg/n_spins_, E_sigma/n_spins_, C, M_avg/n_spins_, M_sigma/n_spins_, Chi);
            }
        }
    }
    
    if(rank_ == 0) {
        display_spins();
        fclose(fptr);
    }
}




double Ising2D::calc_energy(){
    
    double E = 0.0;
    imat nn, nnn;
    
    for(int i = 0; i < N_; ++i){
        for(int j = 0; j < N_; ++j){
            
            nn = nearest_neighbors(i,j);
            nnn = next_nearest_neighbors(i,j);
            
            for(int k = 0; k < 4; ++k){
                E += K_(0)*spins_(i,j)*spins_(nn(k,0),nn(k,1));
                E += K_(1)*spins_(i,j)*spins_(nnn(k,0),nnn(k,1));
            }
        }
    }
    
    return E;
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
            
            for(int k = 0; k < 4; ++k){
                S(0) += spins_(i,j)*spins_(nn(k,0),nn(k,1));
                S(1) += spins_(i,j)*spins_(nnn(k,0),nnn(k,1));
            }
        }
    }
    
    return S;
}

bool Ising2D::in_cluster(int k){
    
    return std::find(cluster_.begin(), cluster_.end(), k) != cluster_.end();
}

void Ising2D::sample_spins(){

    cluster_.clear();
    bool cluster_grew;
    
    // add spins to cluster until either
    // cluster stops growing or cluster is empty
    do{
        cluster_grew = grow_cluster();
    }while(cluster_grew == true && cluster_.size() > 0);
}


bool Ising2D::grow_cluster(){
    
    bool cluster_grew = false;
    int n_cluster = cluster_.size();
    int i, j, k, n;
    double P1 = 1-exp(2*K_(0));
    double P2 = 1-exp(2*K_(1));
    imat nn, nnn;

    // if cluster is empty, add random lattice site
    if(n_cluster == 0){
        
        k = rand_spin_index_(rng);
        cluster_.push_back(k);
        cluster_grew = true;
    }
    
    // if cluster is not empty
    else{
        
        // loop thru elements of current cluster
        for(int c = 0; c < n_cluster; ++c){
            
            k = cluster_[c];
            i = ((k%N_)+N_)%N_;
            j = floor((double)k/N_);
            
            // look thru neighbors of element
            nn = nearest_neighbors(i, j);
            nnn = next_nearest_neighbors(i, j);
            for(int l = 0; l < 4; ++l){
                
                // if spin matches
                if(spins_(i,j) == spins_(nn(l,0), nn(l,1))){
                    
                    // neighbor spin index
                    n = nn(l,1)*N_+nn(l,0);
                    
                    // if neighbor isn't already in cluster
                    if(!in_cluster(n)){
                        
                        // add neighbor to cluster with probability P1
                        if(rand_unif() < P1){
                            cluster_.push_back(n);
                            cluster_grew = true;
                        }
                    }
                }
                
                /*
                // if spin matches
                if(spins_(i,j) == spins_(nnn(l,0), nnn(l,1))){
                    
                    // neighbor spin index
                    n = nnn(l,1)*N_+nnn(l,0);
                    
                    // if neighbor isn't already in cluster
                    if(!in_cluster(n)){
                        
                        // add neighbor to cluster with probability P1
                        if(rand_unif() < P2){
                            cluster_.push_back(n);
                            cluster_grew = true;
                        }
                    }
                }
                 */
            }
            
            // flip current spin
            spins_(i,j) *= -1;
        }
        
        // remove checked elements
        cluster_.erase(cluster_.begin(), cluster_.begin()+n_cluster);
    }

    // return true if spins were added to the cluster
    return cluster_grew;
}


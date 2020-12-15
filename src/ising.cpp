#include "ising.hpp"

IsingModel::IsingModel(vec2D K){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

    K_ = K;
}

void IsingModel::equilibrate(std::shared_ptr<Lattice> pLattice, int n_samples_eq, bool write){

    FILE* fptr = NULL;
    if(rank_ == 0){
        
        printf("Equilibrating %i N = %i lattice(s) at K = ", n_processes_, pLattice->N_);
        print_vec2D(K_);
        
        if(write){
            std::string filename = "equilibrate_N_"+std::to_string(pLattice->N_)
                                   +"_K1_"+get_rounded_str(K_(0))
                                   +"_K2_"+get_rounded_str(K_(1))
                                   +".txt";

            fptr = fopen(filename.c_str(), "w");
            fprintf(fptr, "# Using %i parallel processes\n", n_processes_);
            fprintf(fptr, "# %s, %s, %s, %s, %s, %s, %s\n",
                    "Iteration", "Avg E/spin", "Stddev E/spin", "Heat Capacity",
                    "Avg |M|/spin", "Stddev |M|/spin", "Susceptibility");
        }
    }

    double E, M, E2, M2;
    double E_avg, M_avg, E2_avg, M2_avg;
    double E_sigma, M_sigma;
    double C, Chi;
    
    for(int n = 0; n <= n_samples_eq; ++n){

        sample_new_configuration(pLattice);
        
        if(write && write_iter(n)){
            
            E = calc_energy(pLattice);
            M = abs(calc_magnetization(pLattice));
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
                
                fprintf(fptr, "%i, %10.7e, %10.7e, %10.7e, %10.7e, %10.7e, %10.7e\n",
                        n, E_avg/pLattice->n_spins_, E_sigma/pLattice->n_spins_, C, M_avg/pLattice->n_spins_, M_sigma/pLattice->n_spins_, Chi);
            }
        }
    }
    
    if(rank_ == 0) {
        pLattice->display_spins();
        if(write) fclose(fptr);
    }
}

double IsingModel::calc_energy(std::shared_ptr<Lattice> pLattice){
    
    double E = 0.0;
    imat nn, nnn;
    
    for(int i = 0; i < pLattice->N_; ++i){
        for(int j = 0; j < pLattice->N_; ++j){
            
            nn = pLattice->nearest_neighbors(i,j);
            nnn = pLattice->next_nearest_neighbors(i,j);
            
            for(int k = 0; k < 4; ++k){
                E += K_(0) * pLattice->spins_(i,j) * pLattice->spins_(nn(k,0),nn(k,1));
                E += K_(1) * pLattice->spins_(i,j) * pLattice->spins_(nnn(k,0),nnn(k,1));
            }
        }
    }
    
    return E;
}


double IsingModel::calc_magnetization(std::shared_ptr<Lattice> pLattice){
    
    return pLattice->spins_.sum();
}


void IsingModel::sample_new_configuration(std::shared_ptr<Lattice> pLattice){

    cluster_.clear();
    bool cluster_grew;
    
    // add spins to cluster until either
    // cluster stops growing or cluster is empty
    do{
        cluster_grew = grow_cluster(pLattice);
    }while(cluster_grew == true && cluster_.size() > 0);
}


bool IsingModel::grow_cluster(std::shared_ptr<Lattice> pLattice){
    
    bool cluster_grew = false;
    int i, j, k, n;
    double P1 = 1-exp(2*K_(0));
    double P2 = 1-exp(2*K_(1));
    imat nn, nnn;

    // if cluster is empty, add random lattice site
    int n_cluster = cluster_.size();
    if(n_cluster == 0){
        
        k = pLattice->choose_random_spin();
        cluster_.push_back(k);
        cluster_grew = true;
    }
    
    // if cluster is not empty
    else{
        
        // loop thru elements of current cluster
        for(int c = 0; c < n_cluster; ++c){
            
            k = cluster_[c];
            i = ((k%pLattice->N_)+pLattice->N_)%pLattice->N_;
            j = floor((double)k/pLattice->N_);
            
            // look thru neighbors of element
            nn = pLattice->nearest_neighbors(i, j);
            nnn = pLattice->next_nearest_neighbors(i, j);
            for(int l = 0; l < 4; ++l){
                
                // if spin matches nearest neighbor
                if(pLattice->spins_(i,j) == pLattice->spins_(nn(l,0), nn(l,1))){
                    
                    // neighbor spin index
                    n = nn(l,1)*pLattice->N_+nn(l,0);
                    
                    // if neighbor isn't already in cluster
                    if(!in_cluster(n)){
                        
                        // add neighbor to cluster with probability P1
                        if(rand_unif() < P1){
                            cluster_.push_back(n);
                            cluster_grew = true;
                        }
                    }
                }
                
                
                // if spin matches next nearest neighbor
                if(pLattice->spins_(i,j) == pLattice->spins_(nnn(l,0), nnn(l,1))){
                    
                    // neighbor spin index
                    n = nnn(l,1)*pLattice->N_+nnn(l,0);
                    
                    // if neighbor isn't already in cluster
                    if(!in_cluster(n)){
                        
                        // add neighbor to cluster with probability P1
                        if(rand_unif() < P2){
                            cluster_.push_back(n);
                            cluster_grew = true;
                        }
                    }
                }
                 
            }
            
            // flip current spin
            pLattice->spins_(i,j) *= -1;
        }
        
        // remove checked elements
        cluster_.erase(cluster_.begin(), cluster_.begin()+n_cluster);
    }

    // return true if spins were added to the cluster
    return cluster_grew;
}

bool IsingModel::in_cluster(int k){
    
    return std::find(cluster_.begin(), cluster_.end(), k) != cluster_.end();
}

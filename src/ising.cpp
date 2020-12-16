#include "ising.hpp"

IsingModel::IsingModel(double K){
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

    K_ = K;
    P_ = 1.0-exp(2.0*K_);
    fptr_ = NULL;
}


void IsingModel::equilibrate(std::shared_ptr<Lattice> pLattice,
                             int n_samples_eq,
                             bool write){

    if(rank_ == 0){
        
        printf("Equilibrating %i N = %i lattice(s) at K = %lf\n",
               n_processes_, pLattice->N_, K_);
        
        if(write){
            std::string filename = "equilibrate_N_"+std::to_string(pLattice->N_)
                                   +"_K_"+get_rounded_str(K_)
                                   +".txt";

            fptr_ = fopen(filename.c_str(), "w");
            fprintf(fptr_, "# Nearest neighbor coupling K = %lf\n", K_);
            fprintf(fptr_, "# Temperature T = %lf\n", -1/K_);
            fprintf(fptr_, "# Number of lattice sites = %i\n", pLattice->N_*pLattice->N_);
            fprintf(fptr_, "# Lattice spacing = %i\n", pLattice->a_);
            fprintf(fptr_, "# Using %i parallel processes\n", n_processes_);
            fprintf(fptr_, "# %s, %s, %s, %s, %s, %s, %s\n",
                    "Iteration", "Avg E/spin", "Stddev E/spin", "Heat Capacity",
                    "Avg |M|/spin", "Stddev |M|/spin", "Susceptibility");
        }
    }

    for(int n = 1; n <= n_samples_eq; ++n){

        sample_new_configuration(pLattice);
        
        if(write && write_iter(n)){
            
            E_ = calc_energy(pLattice);
            M_ = abs(calc_magnetization(pLattice));
            E2_ = E_*E_;
            M2_ = M_*M_;
            
            MPI_Reduce(&E_, &E_avg_, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&M_, &M_avg_, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&E2_, &E2_avg_, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&M2_, &M2_avg_, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if(rank_ == 0){
                
                E_avg_ /= n_processes_;
                M_avg_ /= n_processes_;
                E2_avg_ /= n_processes_;
                M2_avg_ /= n_processes_;
                
                E_sigma_ = E2_avg_-E_avg_*E_avg_;
                M_sigma_ = M2_avg_-M_avg_*M_avg_;
                
                C_ = E_sigma_*K_*K_;
                Chi_ = -M_sigma_*K_;
                
                E_sigma_ = sqrt(E_sigma_);
                M_sigma_ = sqrt(M_sigma_);
                
                fprintf(fptr_, "%i, %10.7e, %10.7e, %10.7e, %10.7e, %10.7e, %10.7e\n",
                        n, E_avg_/pLattice->n_spins_, E_sigma_/pLattice->n_spins_, C_, M_avg_/pLattice->n_spins_, M_sigma_/pLattice->n_spins_, Chi_);
            }
        }
    }
    
    if(rank_ == 0){
        pLattice->display_spins();
        if(write){
            pLattice->write_spins(fptr_);
            fclose(fptr_);
        }
    }
}


void IsingModel::sample_new_configuration(std::shared_ptr<Lattice> pLattice){

    // add spins to cluster until cluster stops growing (becomes empty)
    do{
        grow_cluster(pLattice);
    }while(cluster_grew_ == true && cluster_.size() > 0);
}


void IsingModel::grow_cluster(std::shared_ptr<Lattice> pLattice){
    
    cluster_grew_ = false;
    n_cluster_ = cluster_.size();
    int i, j, k, n;
    
    // if cluster is empty, add random lattice site
    if(n_cluster_ == 0){
        
        k = pLattice->choose_random_spin();
        cluster_.push_back(k);
        cluster_grew_ = true;
    }
    
    // if cluster is not empty
    else{
        
        // loop thru elements of current cluster
        for(int c = 0; c < n_cluster_; ++c){
            
            k = cluster_[c];
            i = ((k%pLattice->N_)+pLattice->N_)%pLattice->N_;
            j = floor((double)k/pLattice->N_);
            
            // look thru neighbors of element
            nn_ = pLattice->nearest_neighbors(i, j);
            for(int l = 0; l < 4; ++l){
                
                // if spin matches nearest neighbor
                if(pLattice->spins_(i,j) == pLattice->spins_(nn_(l,0), nn_(l,1))){
                    
                    // neighbor spin index
                    n = nn_(l,1)*pLattice->N_+nn_(l,0);
                    
                    // if neighbor isn't already in cluster
                    if(!in_cluster(n)){
                        
                        // add neighbor to cluster with probability P1
                        if(rand_unif() < P_){
                            cluster_.push_back(n);
                            cluster_grew_ = true;
                        }
                    }
                }
            }
            
            // flip current spin
            pLattice->spins_(i,j) *= -1;
        }
        
        // remove checked cluster elements
        cluster_.erase(cluster_.begin(), cluster_.begin()+n_cluster_);
    }
}


bool IsingModel::in_cluster(int k){
    
    return std::find(cluster_.begin(), cluster_.end(), k) != cluster_.end();
}


double IsingModel::calc_energy(std::shared_ptr<Lattice> pLattice){
    
    E_ = 0.0;
    for(int i = 0; i < pLattice->N_; ++i){
        for(int j = 0; j < pLattice->N_; ++j){
            
            nn_ = pLattice->nearest_neighbors(i,j);
            
            for(int k = 0; k < 4; ++k){
                E_ += K_ * pLattice->spins_(i,j) * pLattice->spins_(nn_(k,0),nn_(k,1));
            }
        }
    }
    
    return E_;
}


double IsingModel::calc_magnetization(std::shared_ptr<Lattice> pLattice){
    
    return pLattice->spins_.sum();
}



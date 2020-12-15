#include "mcrg.hpp"

MonteCarloRenormalizationGroup::MonteCarloRenormalizationGroup(int b){
    
    b_ = b;
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    if(rank_ == 0){
        
        print_bold("\n=============================================\n");
        print_bold("==========       MONTE CARLO       ==========\n");
        print_bold("==========  RENORMALIZATION GROUP  ==========\n");
        print_bold("=============================================\n\n");
        print_bold("* Using "+std::to_string(n_processes_)+" parallel process(es)\n");
        print_bold("* Scaling factor b = "+std::to_string(b_)+"\n");
    }
}

/*
void MonteCarloRenormalizationGroup::calc_critical_exponent(int n_samples_eq,
                                                            int n_samples,
                                                            int N0,
                                                            vec2D Kc){
    // open file
    FILE* fptr = NULL;
    if(rank_ == 0){
        
        print_bold("* Calculating critical exponent at Kc = ");
        print_vec2D(Kc);
        
        
        std::string filename = "critical_exponent_N_"+std::to_string(N0)
                               +"_K1_"+get_rounded_str(Kc(0))
                               +"_K2_"+get_rounded_str(Kc(1))+".txt";
        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# %23s  %25s\n", "Blocking Level n", "Critical Exponent nu");
    }
    
    int n_samples_loc = split_samples(n_samples);
    int n_transformations = floor(log(N0)/log(b_))-1;
    vec2D S, Sb;
    mat S_avg(n_transformations, 2);

    
    // initialize Ising model at critical coupling
    pIsing = new IsingModel(Kc);

    // equilibrate initial lattice at critical coupling
    pLattice = new Lattice(N0);
    pIsing->equilibrate(pLattice, n_samples_eq, false);
    
    // start sampling
    if(rank_ == 0) printf("Sampling %i configurations...\n", n_samples);
    for(int samples = 0; samples < n_samples_loc; ++samples){
        
        // get new configuration
        pIsing->sample_new_configuration(pLattice);
        
        // calculate spin interactions
        S = pLattice->calc_interactions();
        
        //
        S_avg.row(0) =
        
        for(int n = 1; n <= n_transformations; ++n){
            
            // apply block spin transformations on new configuration
            pLatticeb = block_spin_transformation(pLattice);
            
        }
    }

    
    
    if(rank_ == 0){

        fclose(fptr);
    }
}
*/

vec2D MonteCarloRenormalizationGroup::locate_critical_point(int n_iterations,
                                                            int n_samples_eq,
                                                            int n_samples,
                                                            int L,
                                                            vec2D K0){
    FILE* fptr = NULL;
    if(rank_ == 0){
        
        print_bold("* Locating critical point starting from K0 = ");
        print_vec2D(K0);
    
        // open file
        std::string filename = "critical_point_L_"+std::to_string(L)
                               +"_K1_"+get_rounded_str(K0(0))+".txt";
        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# %13s  %15s\n", "Iteration", "K1");
    }
    
    // iteratively approach critical point
    vec2D K = K0;
    for(int i = 1; i <= n_iterations; ++i){
        
        // print iteration
        if(rank_ == 0){
            printf("\nIteration %i: \n", i);
            fprintf(fptr, "%15i, %15.10lf\n", i, K(0));
        }

        // get new approximation
        K = approx_critical_point(n_samples_eq, n_samples, L, K, fptr);
    }
    
    // print results
    if(rank_ == 0){

        fclose(fptr);
        
        print_bold("* Critical point: Kc = ");
        print_vec2D(K);
        
        print_bold("* Critical temperature: Tc = ");
        printf("%.5lf\n", -1.0/K(0));
    }

    return K;
}


vec2D MonteCarloRenormalizationGroup::approx_critical_point(int n_samples_eq,
                                                            int n_samples,
                                                            int L,
                                                            vec2D K,
                                                            FILE* fptr){
    int n_samples_loc = split_samples(n_samples);
    int n_transformations = floor(log(L)/log(b_))-1;
    
    // initialize Ising model at K
    IsingModel* pIsing = new IsingModel(K);
    
    // equilibrate large lattice
    Lattice* pLatticeL = new Lattice(L);
    pIsing->equilibrate(pLatticeL, n_samples_eq, false);
    
    // equilibrate small lattice
    int S = L/b_;
    Lattice* pLatticeS = new Lattice(S);
    pIsing->equilibrate(pLatticeS, n_samples_eq, false);
    
    // containers
    double SL, SL_avg, SL_avg_loc;
    double SS, SS_avg, SS_avg_loc;
    double SLb, SSb;
    vec SLb_avg(n_transformations);
    vec SSb_avg(n_transformations);
    vec SLb_SL_avg(n_transformations);
    vec SSb_SS_avg(n_transformations);
    vec SLb_avg_loc(n_transformations);
    vec SSb_avg_loc(n_transformations);
    vec SLb_SL_avg_loc(n_transformations);
    vec SSb_SS_avg_loc(n_transformations);
    Lattice* pLatticeLb = NULL;
    Lattice* pLatticeSb = NULL;
    
    // set local sums to zero
    SL_avg_loc = 0.0;
    SS_avg_loc = 0.0;
    SLb_avg_loc.setZero();
    SSb_avg_loc.setZero();
    SLb_SL_avg_loc.setZero();
    SSb_SS_avg_loc.setZero();
    
    if(rank_ == 0) printf("Sampling...\n");
    for(int samples = 0; samples < n_samples_loc; ++samples){
        
        // get new configurations
        pIsing->sample_new_configuration(pLatticeL);
        pIsing->sample_new_configuration(pLatticeS);
        
        // calculate nearest neighbor interactions
        SL = pLatticeL->calc_nearest_neighbor_interaction();
        SS = pLatticeS->calc_nearest_neighbor_interaction();
        
        // add up values for averages
        SL_avg_loc += SL;
        SS_avg_loc += SS;
        
        // apply 1 transformation to larger lattice
        pLatticeLb = block_spin_transformation(pLatticeL);
        pLatticeSb = pLatticeS;
        
        // apply n transformations to both lattices
        for(int n = 0; n < n_transformations; ++n){
            
            pLatticeLb = block_spin_transformation(pLatticeLb);
            pLatticeSb = block_spin_transformation(pLatticeSb);
            
            // calculate nearest neighbor interactions
            SLb = pLatticeLb->calc_nearest_neighbor_interaction();
            SSb = pLatticeSb->calc_nearest_neighbor_interaction();
            
            // add up values for averages
            SLb_avg_loc(n) += SLb;
            SSb_avg_loc(n) += SSb;
            SLb_SL_avg_loc(n) += SLb*SL;
            SSb_SS_avg_loc(n) += SSb*SS;
        }
    }
    
    MPI_Allreduce(&SL_avg_loc, &SL_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&SS_avg_loc, &SS_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SLb_avg_loc.data(), SLb_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SSb_avg_loc.data(), SSb_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SLb_SL_avg_loc.data(), SLb_SL_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SSb_SS_avg_loc.data(), SSb_SS_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    SL_avg /= n_samples;
    SS_avg /= n_samples;
    SLb_avg /= n_samples;
    SSb_avg /= n_samples;
    SLb_SL_avg /= n_samples;
    SSb_SS_avg /= n_samples;
    
    // calculate distance to critical point for various blocking levels
    vec2D Kc;
    double dK1, dSL_dK, dSS_dK;
    for(int n = 0; n < n_transformations; ++n){
        
        dSL_dK = SLb_SL_avg(n) - SLb_avg(n) * SL_avg;
        dSS_dK = SSb_SS_avg(n) - SSb_avg(n) * SS_avg;
        dK1 = (SLb_avg(n) - SSb_avg(n)) / (dSL_dK - dSS_dK);
        Kc = K;
        Kc(0) -= dK1;
        
        if(rank_ == 0){
            
            printf("n = %i: Kc = ", n);
            print_vec2D(Kc);
            printf("\n");
        }
    }
    
    // return result for largest n
    return Kc;
}



Lattice* MonteCarloRenormalizationGroup::block_spin_transformation(Lattice* pLattice){
    
    int Nb = pLattice->N_/b_;
    int spin_tot;
    imat block_spins(Nb, Nb);
    
    // loop thru blocks
    for(int ib = 0; ib < Nb; ++ib){
        for(int jb = 0; jb < Nb; ++jb){
            
            // loop thru spins in each block
            spin_tot = 0;
            for(int i = ib*b_; i < (ib+1)*b_; ++i){
                for(int j = jb*b_; j < (jb+1)*b_; ++j){
                    spin_tot += pLattice->spins_(i,j);
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
    
    Lattice* pLatticeb = new Lattice(pLattice->a_*b_, block_spins);

    return pLatticeb;
}



int MonteCarloRenormalizationGroup::split_samples(int n_samples){
    
    int n_samples_loc = ceil((double)n_samples/(double)n_processes_);
    if(rank_ == 0){
        n_samples_loc = n_samples-(n_processes_-1)*n_samples_loc;
    }
    
    return n_samples_loc;
}

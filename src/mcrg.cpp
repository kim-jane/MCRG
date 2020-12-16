#include "mcrg.hpp"

MonteCarloRenormalizationGroup::MonteCarloRenormalizationGroup(int b){
    
    b_ = b;
    fptr_ = NULL;
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    if(rank_ == 0){
        
        printf("\n=============================================\n");
        printf("==========       MONTE CARLO       ==========\n");
        printf("==========  RENORMALIZATION GROUP  ==========\n");
        printf("=============================================\n\n");
        printf("* Using %i parallel process(es)\n", n_processes_);
        printf("* Scaling factor b = %i\n", b_);
    }
}


void MonteCarloRenormalizationGroup::calc_critical_exponent(int n_samples_eq,
                                                            int n_samples,
                                                            int N,
                                                            double K){
    // open file
    if(rank_ == 0){
        
        printf("* Calculating critical exponent at K = %lf\n\n", K);
        
        std::string filename = "critical_exponent_N_"+std::to_string(N)
                               +"_K_"+get_rounded_str(K)+".txt";
        fptr_ = fopen(filename.c_str(), "w");
        
        fprintf(fptr_, "# Number of parallel processes = %i\n", n_processes_);
        fprintf(fptr_, "# Number of equilibration samples = %i\n", n_samples_eq);
        fprintf(fptr_, "# Number of samples = %i\n", n_samples);
        fprintf(fptr_, "# %23s  %25s  %25s\n",
                "Blocking Level n", "Largest Eigenvalue", "Critical Exponent nu");
    }
    
    int n_samples_loc = split_samples(n_samples);
    int n_transformations = floor(log(N)/log(b_));

    // initialize Ising model at K
    std::unique_ptr<IsingModel> pIsing(new IsingModel(K));
    
    // equilibrate lattice
    std::shared_ptr<Lattice> pLattice(new Lattice(N));
    pIsing->equilibrate(pLattice, n_samples_eq, false);
    
    // containers
    vec S_vec, Sb_vec;
    vec Sb_S_flattened, Sb_Sb_flattened;
    mat Sb_S, Sb_Sb;
    mat S(n_transformations+1,2);
    mat S_avg(n_transformations+1,2);
    mat Sb_S_avg(n_transformations,4);
    mat Sb_Sb_avg(n_transformations,4);
    mat S_avg_loc(n_transformations+1,2);
    mat Sb_S_avg_loc(n_transformations,4);
    mat Sb_Sb_avg_loc(n_transformations,4);
    std::shared_ptr<Lattice> pLatticeb;
    
    // set local sums to zero
    S_avg_loc.setZero();
    Sb_S_avg_loc.setZero();
    Sb_Sb_avg_loc.setZero();
    
    // start sampling
    if(rank_ == 0) printf("Sampling %i configurations...\n", n_samples);
    for(int samples = 0; samples < n_samples_loc; ++samples){
        
        // get new configuration
        pIsing->sample_new_configuration(pLattice);
        
        // calculate spin interactions at each blocking level
        S.row(0) = pLattice->calc_interactions();
        pLatticeb = pLattice;
        for(int n = 1; n <= n_transformations; ++n){
            
            pLatticeb = block_spin_transformation(pLatticeb);
            S.row(n) = pLatticeb->calc_interactions();
            
            // add up values for averages
            Sb_vec = S.row(n);
            S_vec = S.row(n-1);
            Sb_S = Sb_vec*S_vec.transpose();
            Sb_Sb = Sb_vec*Sb_vec.transpose();
            Sb_S_flattened = flatten(Sb_S);
            Sb_Sb_flattened = flatten(Sb_Sb);
            Sb_S_avg_loc.row(n-1) += Sb_S_flattened;
            Sb_Sb_avg_loc.row(n-1) += Sb_Sb_flattened;
        }
        
        // add up values for averages
        S_avg_loc += S;
    }
    
    // sum results from all parallel processes
    MPI_Allreduce(S_avg_loc.data(), S_avg.data(), 2*(n_transformations+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(Sb_S_avg_loc.data(), Sb_S_avg.data(), 4*n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(Sb_Sb_avg_loc.data(), Sb_Sb_avg.data(), 4*n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // average
    S_avg /= n_samples;
    Sb_S_avg /= n_samples;
    Sb_Sb_avg /= n_samples;
    
    // calculate RG eigenvalues at various blocking levels
    mat2D dSb_dK, dSb_dKb, T;
    double lambda, lambda1, lambda2, nu;
    for(int n = 0; n < n_transformations; ++n){
        
        Sb_S_flattened = Sb_S_avg.row(n);
        Sb_Sb_flattened = Sb_Sb_avg.row(n);
        Sb_vec = S_avg.row(n+1);
        S_vec = S_avg.row(n);
        
        dSb_dK = unflatten(Sb_S_flattened)-Sb_vec*S_vec.transpose();
        dSb_dKb = unflatten(Sb_Sb_flattened)-Sb_vec*Sb_vec.transpose();
        T = dSb_dK * dSb_dKb.inverse();
        
        // get largest eigenvalue
        EigenSolver<mat2D> solver(T);
        lambda1 = solver.eigenvalues()(0).real();
        lambda2 = solver.eigenvalues()(1).real();
        if(lambda1 > lambda2) lambda = lambda1;
        else lambda = lambda2;
        nu = log(lambda)/log(b_);
        
        if(rank_ == 0){
            printf("n = %i: lambda = %lf, nu = %lf\n", n, lambda, nu);
            fprintf(fptr_, "%25i, %25.10lf, %25.10lf\n", n, lambda, nu);
        }
    }
    
    if(rank_ == 0){
        fclose(fptr_);
    }
}

double MonteCarloRenormalizationGroup::locate_critical_point(int n_iterations,
                                                             int n_samples_eq,
                                                             int n_samples,
                                                             int L,
                                                             double K0){
    if(rank_ == 0){
        
        // print starting point
        printf("* Locating critical point starting from K0 = %lf\n", K0);
    
        // open file
        std::string filename = "critical_point_L_"+std::to_string(L)
                               +"_K_"+get_rounded_str(K0)+".txt";
        fptr_ = fopen(filename.c_str(), "w");
        
        // headings
        fprintf(fptr_, "# Number of parallel processes = %i\n", n_processes_);
        fprintf(fptr_, "# Number of equilibration samples = %i\n", n_samples_eq);
        fprintf(fptr_, "# Number of samples = %i\n", n_samples);
        fprintf(fptr_, "# %13s  %15s  %15s  %15s\n",
                "Iteration", "Blocking Level n", "Starting K","Approximate Kc");
    }
    
    // iteratively approach critical point
    double K = K0;
    for(iter_ = 1; iter_ <= n_iterations; ++iter_){
        
        if(rank_ == 0) printf("\nIteration %i:\n", iter_);
        
        // get new approximation
        K = approx_critical_point(n_samples_eq, n_samples, L, K);
    }
    
    // print results
    if(rank_ == 0){

        printf("* Critical point: Kc = %lf\n", K);
        printf("* Critical temperature: Tc = %lf\n", -1.0/K);
        fclose(fptr_);
    }
    
    return K;
}


double MonteCarloRenormalizationGroup::approx_critical_point(int n_samples_eq,
                                                             int n_samples,
                                                             int L,
                                                             double K){
    int n_samples_loc = split_samples(n_samples);
    int n_transformations = floor(log(L)/log(b_))-2;
    
    // initialize Ising model at K
    std::unique_ptr<IsingModel> pIsing(new IsingModel(K));
    
    // equilibrate large lattice
    std::shared_ptr<Lattice> pLatticeL(new Lattice(L));
    pIsing->equilibrate(pLatticeL, n_samples_eq, false);
    
    // equilibrate small lattice
    int S = L/b_;
    std::shared_ptr<Lattice> pLatticeS(new Lattice(S));
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
    std::shared_ptr<Lattice> pLatticeLb;
    std::shared_ptr<Lattice> pLatticeSb;
    
    // set local sums to zero
    SL_avg_loc = 0.0;
    SS_avg_loc = 0.0;
    SLb_avg_loc.setZero();
    SSb_avg_loc.setZero();
    SLb_SL_avg_loc.setZero();
    SSb_SS_avg_loc.setZero();
    
    if(rank_ == 0) printf("Sampling %i configurations each...\n", n_samples);
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
            
            // transform
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
    
    // sum results from all parallel processes
    MPI_Allreduce(&SL_avg_loc, &SL_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&SS_avg_loc, &SS_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SLb_avg_loc.data(), SLb_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SSb_avg_loc.data(), SSb_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SLb_SL_avg_loc.data(), SLb_SL_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SSb_SS_avg_loc.data(), SSb_SS_avg.data(), n_transformations, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // average
    SL_avg /= n_samples;
    SS_avg /= n_samples;
    SLb_avg /= n_samples;
    SSb_avg /= n_samples;
    SLb_SL_avg /= n_samples;
    SSb_SS_avg /= n_samples;
    
    // calculate distance to critical point for various blocking levels
    double Kc = 0;
    double dK, dSL_dK, dSS_dK;
    for(int n = 0; n < n_transformations; ++n){
        
        dSL_dK = SLb_SL_avg(n) - SLb_avg(n) * SL_avg;
        dSS_dK = SSb_SS_avg(n) - SSb_avg(n) * SS_avg;
        dK = (SLb_avg(n) - SSb_avg(n)) / (dSL_dK - dSS_dK);
        Kc = K-dK;
        
        if(rank_ == 0){
            
            printf("n = %i: Kc = %lf\n", n, Kc);
            fprintf(fptr_, "%15i, %15i, %15.10lf, %15.10lf\n", iter_, n, K, Kc);
        }
    }
    
    // return result for largest n
    return Kc;
}



std::shared_ptr<Lattice> MonteCarloRenormalizationGroup::block_spin_transformation(std::shared_ptr<Lattice> pLattice){
    
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
    
    std::shared_ptr<Lattice> pLatticeb(new Lattice(pLattice->a_*b_, block_spins));

    return pLatticeb;
}



int MonteCarloRenormalizationGroup::split_samples(int n_samples){
    
    int n_samples_loc = ceil((double)n_samples/(double)n_processes_);
    if(rank_ == 0){
        n_samples_loc = n_samples-(n_processes_-1)*n_samples_loc;
    }
    
    return n_samples_loc;
}

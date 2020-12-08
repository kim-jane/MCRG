#include "mcrg.hpp"

MonteCarloRenormalizationGroup::MonteCarloRenormalizationGroup(int b,
                                                               int verbose){
    
    b_ = b;
    verbose_ = verbose;
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    if(rank_ == 0){
        printf("\n=============================================\n");
        printf("==========       MONTE CARLO       ==========\n");
        printf("==========  RENORMALIZATION GROUP  ==========\n");
        printf("=============================================\n\n");
        
        print_bold("* Using "+std::to_string(n_processes_)+" parallel process(es)\n");
        print_bold("* Scaling factor b = "+std::to_string(b_)+"\n");
        
        if(verbose_ > 1 && n_processes_ > 1){
            print_bold("* Printed spin lattices belong to root process\n");
        }
    }
}

double MonteCarloRenormalizationGroup::calc_critical_exponent(int n_samples,
                                                              int N0,
                                                              vec2D Kc){
    FILE* fptr = NULL;
    if(rank_ == 0){
        
        // print starting point
        print_bold("* Calculating critical exponent nu at Kc = ");
        print_vec2D(Kc);
        
        // open file
        std::string Kstr = get_string(Kc);
        std::string filename = "critical_exponent_s_"+std::to_string(n_samples)
                               +"_N_"+std::to_string(N0)+Kstr+".txt";
        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# %23s  %25s\n", "Blocking Level n", "Critical Exponent nu");
    }
    
    int n_samples_eq = 1E4;
    int n_samples_loc = split_samples(n_samples);
    int n_transformations = floor(log(N0)/log(b_))-2;
    
}


vec2D MonteCarloRenormalizationGroup::locate_critical_point(int n_iterations,
                                                            int n_samples,
                                                            int L,
                                                            int S,
                                                            vec2D K0){
    FILE* fptr = NULL;
    if(rank_ == 0){
        
        // print starting point
        print_bold("* Locating critical point Kc starting from K0 = ");
        print_vec2D(K0);
        
        // open file
        std::string Kstr = get_string(K0);
        std::string filename = "critical_point_s_"+std::to_string(n_samples)
                               +"_L_"+std::to_string(L)
                               +"_S_"+std::to_string(S)+Kstr+".txt";
        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# %13s  %15s  %15s\n", "Iteration", "K1", "K2");
    }
    
    
    vec2D K = K0;
    for(int i = 1; i <= n_iterations; ++i){
        
        if(rank_ == 0){
            
            printf("\nIteration %i:\n", i);
            
            fprintf(fptr, "%15i, %15.10lf, %15.10lf,\n", i, K(0), K(1));
            fclose(fptr);
        }

        // get new estimate of critical point
        K = approx_critical_point(n_samples, L, S, K);
    }

    if(rank_ == 0){
        
        fprintf(fptr, "%15i, %15.10lf, %15.10lf,\n", n_iterations, K(0), K(1));
        fclose(fptr);
        
        print_bold("* Critical point: Kc = ");
        print_vec2D(K);
        
        print_bold("* Critical temperature: Tc = ");
        printf("%.5lf\n", -1.0/K(0));
    }

    return K;
}


vec2D MonteCarloRenormalizationGroup::approx_critical_point(int n_samples,
                                                            int L,
                                                            int S,
                                                            vec2D K){
    
    FILE* fptr = NULL;
    int n_samples_eq = 1E4;
    int n_samples_loc = split_samples(n_samples);
    int n_transformations = floor(log(S)/log(b_));
    int m = floor(log(L/S)/log(b_));

    // equilibrate initial large lattice
    Ising2D* pIsingL;
    pIsingL = new Ising2D(L, K, verbose_);
    pIsingL->equilibrate(n_samples_eq, fptr);
    
    // apply m transformations
    Ising2D* pIsingLm = pIsingL;
    for(int n = 0; n < m; ++n){
        pIsingLm = pIsingLm->block_spin_transformation(b_);
    }
    
    // equilibrate small lattice with the same number of
    // lattice sites as transformed large lattice
    Ising2D* pIsingS = new Ising2D(pIsingLm->N_, K, verbose_);
    pIsingS->equilibrate(n_samples_eq, fptr);
    
    // transformed small lattice at n = 0
    Ising2D* pIsingSm = pIsingS;
    
    vec2D Kc, dK;
    vec2D SL, SLm, SS, SSm;
    vec2D SL_avg, SLm_avg, SS_avg, SSm_avg;
    vec2D SL_avg_loc, SLm_avg_loc, SS_avg_loc, SSm_avg_loc;
    mat2D SLm_SL_avg, SSm_SS_avg;
    mat2D SLm_SL_avg_loc, SSm_SS_avg_loc;
    
    for(int n = 0; n < n_transformations; ++n){
    
        // calculate correlation functions
        SL_avg_loc.setZero();
        SLm_avg_loc.setZero();
        SS_avg_loc.setZero();
        SSm_avg_loc.setZero();
        SLm_SL_avg_loc.setZero();
        SSm_SS_avg_loc.setZero();
        
        for(int samples = 0; samples < n_samples_loc; ++samples){
            
            // get new samples
            pIsingL->sample_spins();
            pIsingLm->sample_spins();
            pIsingS->sample_spins();
            pIsingSm->sample_spins();
            
            // calculate NN and NNN interactions for each lattice
            SL = pIsingL->calc_spin_interactions();
            SLm = pIsingLm->calc_spin_interactions();
            SS = pIsingS->calc_spin_interactions();
            SSm = pIsingSm->calc_spin_interactions();
            
            // add up values for averages
            SL_avg_loc += SL;
            SLm_avg_loc += SLm;
            SS_avg_loc += SS;
            SSm_avg_loc += SSm;
            SLm_SL_avg_loc += SLm*SL.transpose();
            SSm_SS_avg_loc += SSm*SS.transpose();
        }
        
        MPI_Allreduce(SL_avg_loc.data(), SL_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SLm_avg_loc.data(), SLm_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SS_avg_loc.data(), SS_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SSm_avg_loc.data(), SSm_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SLm_SL_avg_loc.data(), SLm_SL_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SSm_SS_avg_loc.data(), SSm_SS_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        SL_avg /= n_samples;
        SLm_avg /= n_samples;
        SS_avg /= n_samples;
        SSm_avg /= n_samples;
        SLm_SL_avg /= n_samples;
        SSm_SS_avg /= n_samples;
        
        mat2D dSL_dK = SLm_SL_avg-SLm_avg*SL_avg.transpose();
        mat2D dSS_dK = SSm_SS_avg-SSm_avg*SS_avg.transpose();
        vec2D dK = (dSL_dK-dSS_dK).inverse() * (SLm_avg-SSm_avg);
        Kc = K-dK;
        
        if(rank_ == 0 && verbose_ > 0){
            
            printf("Kc = ");
            print_vec2D(Kc);
        }
        if(n < n_transformations-1){
            pIsingLm = pIsingLm->block_spin_transformation(b_);
            pIsingSm = pIsingSm->block_spin_transformation(b_);
        }
    }
    
    return Kc;
}


std::string MonteCarloRenormalizationGroup::get_string(vec2D K){
    
    std::stringstream ss_K1, ss_K2;
    ss_K1 << std::setprecision(3) << K(0);
    ss_K2 << std::setprecision(3) << K(1);
    
    std::string str = "_K1_"+ss_K1.str()
                      +"_K2_"+ss_K2.str();
    
    return str;
}

int MonteCarloRenormalizationGroup::split_samples(int n_samples){
    
    int n_samples_loc = ceil((double)n_samples/(double)n_processes_);
    if(rank_ == 0){
        n_samples_loc = n_samples-(n_processes_-1)*n_samples_loc;
    }
    
    return n_samples_loc;
}

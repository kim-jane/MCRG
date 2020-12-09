#include "mcrg.hpp"

MonteCarloRenormalizationGroup::MonteCarloRenormalizationGroup(int b){
    
    b_ = b;
    
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    
    if(rank_ == 0){
        
        printf("\n=============================================\n");
        printf("==========       MONTE CARLO       ==========\n");
        printf("==========  RENORMALIZATION GROUP  ==========\n");
        printf("=============================================\n\n");
        
        print_bold("* Using "+std::to_string(n_processes_)+" parallel process(es)\n");
        print_bold("* Scaling factor b = "+std::to_string(b_)+"\n");
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
        std::string filename = "critical_exponent_s_"+std::to_string(n_samples)
                               +"_N_.txt";
        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# %23s  %25s\n", "Blocking Level n", "Critical Exponent nu");
    }
    
    int n_samples_eq = 1E5;
    int n_samples_loc = split_samples(n_samples);
    int n_transformations = floor(log(N0)/log(b_))-1;
    
}


vec2D MonteCarloRenormalizationGroup::locate_critical_point(int n_iterations,
                                                            int n_samples,
                                                            int L0,
                                                            vec2D K0){
    FILE* fptr = NULL;
    if(rank_ == 0){
        
        // print starting point
        print_bold("* Locating critical point Kc starting from K0 = ");
        print_vec2D(K0);
        
        // open file

        std::string filename = "critical_point_s_"+std::to_string(n_samples)
                               +"_L_"+std::to_string(L0)+".txt";
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
        K = approx_critical_point(n_samples, L0, K);
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
                                                            int L0,
                                                            vec2D K){
    
    int S0 = L0/b_;
    int n_transformations = floor(log(S0)/log(b_))-2;
    int n_samples_loc = split_samples(n_samples);
    int n_samples_eq = 1E5;

    // equilibrate initial large lattice
    Ising2D* pIsingL0;
    pIsingL0 = new Ising2D(L0, K);
    pIsingL0->equilibrate(n_samples_eq, false);
    
    // apply 1 transformation to large lattice
    Ising2D* pIsingL = pIsingL0->block_spin_transformation(b_);
    
    // equilibrate small lattice with the same number of
    // lattice sites as transformed large lattice
    Ising2D* pIsingS0 = new Ising2D(S0, K);
    pIsingS0->equilibrate(n_samples_eq, false);
    
    // transformed small lattice at n = 0
    Ising2D* pIsingS = pIsingS0;
    
    vec2D Kc, dK;
    vec2D SL0, SS0, SL, SS;
    vec2D SL0_avg, SS0_avg, SL_avg, SS_avg;
    vec2D SL0_avg_loc, SS0_avg_loc, SL_avg_loc, SS_avg_loc;
    mat2D SL_SL0_avg, SS_SS0_avg;
    mat2D SL_SL0_avg_loc, SS_SS0_avg_loc;
    
    for(int n = 0; n < n_transformations; ++n){
    
        // calculate correlation functions
        SL0_avg_loc.setZero();
        SS0_avg_loc.setZero();
        SL_avg_loc.setZero();
        SS_avg_loc.setZero();
        SL_SL0_avg_loc.setZero();
        SS_SS0_avg_loc.setZero();
        
        for(int samples = 0; samples < n_samples_loc; ++samples){
            
            // get new samples
            pIsingL0->sample_spins();
            pIsingS0->sample_spins();
            pIsingL->sample_spins();
            pIsingS->sample_spins();
            
            // calculate spin interactions for each lattice
            SL0 = pIsingL0->calc_spin_interactions();
            SS0 = pIsingS0->calc_spin_interactions();
            SL = pIsingL->calc_spin_interactions();
            SS = pIsingS->calc_spin_interactions();
            
            // add up values for averages
            SL0_avg_loc += SL0;
            SS0_avg_loc += SS0;
            SL_avg_loc += SL;
            SS_avg_loc += SS;
            SL_SL0_avg_loc += SL*SL0.transpose();
            SS_SS0_avg_loc += SS*SS0.transpose();
        }
        
        MPI_Allreduce(SL0_avg_loc.data(), SL0_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SS0_avg_loc.data(), SS0_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SL_avg_loc.data(), SL_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SS_avg_loc.data(), SS_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SL_SL0_avg_loc.data(), SL_SL0_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(SS_SS0_avg_loc.data(), SS_SS0_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        SL0_avg /= n_samples;
        SS0_avg /= n_samples;
        SL_avg /= n_samples;
        SS_avg /= n_samples;
        SL_SL0_avg /= n_samples;
        SS_SS0_avg /= n_samples;
        
        mat2D dSL_dK = SL_SL0_avg-SL_avg*SL0_avg.transpose();
        mat2D dSS_dK = SS_SS0_avg-SS_avg*SS0_avg.transpose();
        vec2D dK = (dSL_dK-dSS_dK).inverse() * (SL_avg-SS_avg);
        Kc = K-dK;
        
        if(rank_ == 0){
            
            printf("Kc = ");
            print_vec2D(Kc);
        }
        if(n < n_transformations-1){
            //pIsingL0 = pIsingL0->block_spin_transformation(b_);
            //pIsingS0 = pIsingS0->block_spin_transformation(b_);
            pIsingL = pIsingL->block_spin_transformation(b_);
            pIsingS = pIsingS->block_spin_transformation(b_);
        }
    }
    
    return Kc;
}




int MonteCarloRenormalizationGroup::split_samples(int n_samples){
    
    int n_samples_loc = ceil((double)n_samples/(double)n_processes_);
    if(rank_ == 0){
        n_samples_loc = n_samples-(n_processes_-1)*n_samples_loc;
    }
    
    return n_samples_loc;
}

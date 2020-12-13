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
    // open file
    FILE* fptr = NULL;
    if(rank_ == 0){
        std::string filename = "critical_exponent_N_"+std::to_string(N0)
                               +"_K1_"+get_rounded_str(Kc(0))
                               +"_K2_"+get_rounded_str(Kc(1))
                               +"_s_"+std::to_string(n_samples)+".txt";
        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# %23s  %25s\n", "Blocking Level n", "Critical Exponent nu");
    }
    
    int n_samples_eq = 1E5;
    int n_samples_loc = split_samples(n_samples);
    
    // get initial estimate of distance to fixed point dK
    //dK = distance_critical_point
    
    
    // equilibrate initial system at critical coupling
    Ising2D* pIsing = new Ising2D(N0, Kc);
    pIsing->equilibrate(n_samples_eq, false);
    if(rank_ == 0) pIsing->display_spins();
    
    // apply one RG transformation
    Ising2D* pIsingb = pIsing->block_spin_transformation(b_);
    
    vec2D S, Sb;
    vec2D S_avg, Sb_avg;
    vec2D S_avg_loc, Sb_avg_loc;
    mat2D Sb_S_avg, Sb_Sb_avg;
    mat2D Sb_S_avg_loc, Sb_Sb_avg_loc;
    mat2D dSb_dK, dSb_dKb;
    mat2D T;
    
    // calculate correlation functions

    S_avg_loc.setZero();
    Sb_avg_loc.setZero();
    Sb_S_avg_loc.setZero();
    Sb_Sb_avg_loc.setZero();
    
    for(int samples = 0; samples < n_samples_loc; ++samples){
        
        // get new samples
        pIsing->sample_spins();
        pIsingb->sample_spins();
        
        // calculate spin interactions for each lattice
        S = pIsing->calc_spin_interactions();
        Sb = pIsingb->calc_spin_interactions();
        
        // add up values for averages
        S_avg_loc += S;
        Sb_avg_loc += Sb;
        Sb_S_avg_loc += Sb*S.transpose();
        Sb_Sb_avg_loc += Sb*Sb.transpose();
    }

    MPI_Allreduce(S_avg_loc.data(), S_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(Sb_avg_loc.data(), Sb_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(Sb_S_avg_loc.data(), Sb_S_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(Sb_Sb_avg_loc.data(), Sb_Sb_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    S_avg /= n_samples;
    Sb_avg /= n_samples;
    Sb_S_avg /= n_samples;
    Sb_Sb_avg /= n_samples;
    
    dSb_dK = Sb_S_avg-Sb_avg*S_avg.transpose();
    dSb_dKb = Sb_Sb_avg-Sb_avg*Sb_avg.transpose();
    T = dSb_dKb.inverse() * dSb_dK;
    EigenSolver<mat2D> solver(T);
    
    std::cout << solver.eigenvalues() << std::endl;

    
    return 0.0;
}


vec2D MonteCarloRenormalizationGroup::locate_critical_point(int n_iterations,
                                                            int n_samples,
                                                            int L,
                                                            vec2D K0){
    FILE* fptr = NULL;
    if(rank_ == 0){
    
        // open file
        std::string filename = "critical_point_L_"+std::to_string(L)
                               +"_K1_"+get_rounded_str(K0(0))
                               +"_K2_"+get_rounded_str(K0(1))
                               +"_s_"+std::to_string(n_samples)+".txt";
        fptr = fopen(filename.c_str(), "w");
        fprintf(fptr, "# %13s  %15s  %15s  %15s  %15s\n", "Iteration", "K1", "K2", "K1'", "K2'");
    }
    
    // iteratively approach critical point
    vec2D K = K0;
    for(int i = 1; i <= n_iterations; ++i){
        
        K(1) = 0.0;
        
        if(rank_ == 0){
            printf("\nIteration %i: Kc = ", i);
            print_vec2D(K);
            fprintf(fptr, "%15i, %15.10lf, %15.10lf, ", i, K(0), K(1));
        }

        K = approx_critical_point(n_samples, L, K);
        
        if(rank_ == 0){
            fprintf(fptr, "%15.10lf, %15.10lf,\n", K(0), K(1));
        }
    }
    
    // print results
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
                                                            vec2D K){
    
    int n_samples_loc = split_samples(n_samples);
    int n_samples_eq = 1E4;
    
    // equilibrate initial large lattice
    Ising2D* pIsingL;
    pIsingL = new Ising2D(L, K);
    pIsingL->equilibrate(n_samples_eq, true);
    
    // apply 1 transformation to large lattice
    Ising2D* pIsingLb = pIsingL->block_spin_transformation(b_);
    
    // equilibrate small lattice with the same number of
    // lattice sites as transformed large lattice
    Ising2D* pIsingS = new Ising2D(L/b_, K);
    pIsingS->equilibrate(n_samples_eq, true);
    
    // containers
    vec2D SL, SL_avg, SL_avg_loc;
    vec2D SLb, SLb_avg, SLb_avg_loc;
    vec2D SS, SS_avg, SS_avg_loc;
    mat2D SLb_SL_avg, SLb_SL_avg_loc;
    mat2D SS_SS_avg, SS_SS_avg_loc;
    
    // calculate correlation functions
    SL_avg_loc.setZero();
    SLb_avg_loc.setZero();
    SS_avg_loc.setZero();
    SLb_SL_avg_loc.setZero();
    SS_SS_avg_loc.setZero();
    
    for(int samples = 0; samples < n_samples_loc; ++samples){
        
        // get new samples
        pIsingL->sample_spins();
        pIsingLb->sample_spins();
        pIsingS->sample_spins();
        
        // calculate spin interactions for each lattice
        SL = pIsingL->calc_spin_interactions();
        SLb = pIsingLb->calc_spin_interactions();
        SS = pIsingS->calc_spin_interactions();
        
        // add up values for averages
        SL_avg_loc += SL;
        SLb_avg_loc += SLb;
        SS_avg_loc += SS;
        SLb_SL_avg_loc += SLb*SL.transpose();
        SS_SS_avg_loc += SS*SS.transpose();
    }
    
    MPI_Allreduce(SL_avg_loc.data(), SL_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SLb_avg_loc.data(), SLb_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SS_avg_loc.data(), SS_avg.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SLb_SL_avg_loc.data(), SLb_SL_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(SS_SS_avg_loc.data(), SS_SS_avg.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    SL_avg /= n_samples;
    SLb_avg /= n_samples;
    SS_avg /= n_samples;
    SLb_SL_avg /= n_samples;
    SS_SS_avg /= n_samples;
    
    mat2D dSL_dK = SLb_SL_avg - SLb_avg * SL_avg.transpose();
    mat2D dSS_dK = SS_SS_avg - SS_avg * SS_avg.transpose();
    vec2D dK = (dSL_dK-dSS_dK).inverse() * (SL_avg-SS_avg);
    vec Kc = K-dK;
    
    return Kc;
}


int MonteCarloRenormalizationGroup::split_samples(int n_samples){
    
    int n_samples_loc = ceil((double)n_samples/(double)n_processes_);
    if(rank_ == 0){
        n_samples_loc = n_samples-(n_processes_-1)*n_samples_loc;
    }
    
    return n_samples_loc;
}

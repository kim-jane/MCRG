#include "mcrg.hpp"

MonteCarloRenormalizationGroup::MonteCarloRenormalizationGroup(int b){
    
    b_ = b;
    verbose_ = true;
    
    printf("\n%10s===================================\n", "");
    printf("%10s=====       MONTE CARLO       =====\n", "");
    printf("%10s=====  RENORMALIZATION GROUP  =====\n", "");
    printf("%10s===================================\n\n", "");
    
    print_bold("* Scaling factor b = "+std::to_string(b)+"\n\n");
    print_bold("* Periodic boundary conditions\n\n");
}

void MonteCarloRenormalizationGroup::run(int n_samples,
                                         int N,
                                         vec2D K){
    
    vec2D K_new = approx_critical_point(n_samples, N, K);
    
    for(int i = 0; i < 2; ++i){
        K_new =  approx_critical_point(n_samples, N, K_new);
    }
    
}

vec2D MonteCarloRenormalizationGroup::approx_critical_point(int n_samples,
                                                            int N,
                                                            vec2D K){
    
    if(N%b_ != 0){
        printf("Lattice size is not divisible by the scaling factor.\n");
        exit(1);
    }
    
    print_bold("* Approximating critical point:\n\n");
    
    // two lattices with sizes differing by b
    Ising2D* pIsing1 = new Ising2D(N, K);
    Ising2D* pIsing2 = new Ising2D(N/b_, K);
    
    // equilibrate both lattices
    pIsing1->equilibrate();
    pIsing2->equilibrate();
    
    // transform larger lattice
    Ising2D* pIsing1b = pIsing1->block_spin_transformation(b_);
    
    // calculate correlation functions
    vec2D S1, S2, S1b;
    vec2D S1_avg, S2_avg, S1b_avg;
    mat2D S1b_S1_avg, S2_S2_avg;
    S1_avg.setZero();
    S2_avg.setZero();
    S1b_avg.setZero();
    S1b_S1_avg.setZero();
    S2_S2_avg.setZero();
    
    printf("Sampling %i configurations each...\n\n", n_samples);
    
    for(int samples = 0; samples < n_samples; ++samples){
        
        // get new samples
        pIsing1->sample_spins();
        pIsing2->sample_spins();
        pIsing1b->sample_spins();
        
        // calculate correlation for each lattice
        S1 = pIsing1->calc_correlation();
        S2 = pIsing2->calc_correlation();
        S1b = pIsing1b->calc_correlation();
        
        // add up values for averages
        S1_avg += S1;
        S2_avg += S2;
        S1b_avg += S1b;
        S1b_S1_avg += S1b*S1.transpose();
        S2_S2_avg += S2*S2.transpose();
    }

    S1_avg /= n_samples;
    S2_avg /= n_samples;
    S1b_avg /= n_samples;
    S1b_S1_avg /= n_samples;
    S2_S2_avg /= n_samples;
    
    mat2D dS1_dK = S1b_S1_avg-S1b_avg*S1_avg.transpose();
    mat2D dS2_dK = S2_S2_avg-S2_avg*S2_avg.transpose();
    
    vec2D dK = (dS1_dK-dS2_dK).inverse() * (S1b_avg-S2_avg);
    vec2D Kc = K-dK;
    
    printf("Distance to critical point: dK = (%.10lf, %.10lf)\n\n", dK(0), dK(1));
    printf("Approximate critical point: Kc = (%.10lf, %.10lf)\n\n", Kc(0), Kc(1));
    
    return Kc;
}

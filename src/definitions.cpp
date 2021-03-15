#include "definitions.hpp"

std::random_device rd;
std::mt19937_64 rng(rd());
std::uniform_int_distribution<int> binary(0, 1);
std::uniform_real_distribution<double> unif(0.0, 1.0);
std::normal_distribution<double> narrow_norm(0.0, 0.1);

vec flatten(mat M){
    
    int n = M.cols();
    vec v(n*n);
    
    for(int j = 0; j < n; ++j){
        v.segment(j*n, n) = M.col(j);
    }
    
    return v;
}

mat unflatten(vec v){
    
    int n = floor(sqrt(v.size()));
    mat M(n,n);
    
    for(int j = 0; j < n; ++j){
        M.col(j) = v.segment(j*n, n);
    }
    
    return M;
}


void display_spin_up(){
    
    printf("\033[34m%2s\033[0m", "O");
}

void display_spin_down(){
    
    printf("\033[37m%2s\033[0m", "X");
}

bool write_iter(int i){
    
    bool write = false;
    
    if(i < 10){
        write = true;
    }
    else if(i < 100 && i%10 == 0){
        write = true;
    }
    else if(i < 1000 && i%100 == 0){
        write = true;
    }
    else if(i < 10000 && i%1000 == 0){
        write = true;
    }
    else if(i < 100000 && i%10000 == 0){
        write = true;
    }
    else if(i%100000 == 0){
        write = true;
    }
    
    return write;
}

std::string get_rounded_str(double num, int precision){
    
    std::stringstream ss;
    ss << std::setprecision(precision) << num;
    std::string str = ss.str();
    
    return str;
}

int split_samples(int rank, int n_processes, int n_samples){
    
    int n_samples_loc = ceil((double)n_samples/(double)n_processes);
    if(rank == 0){
        n_samples_loc = n_samples-(n_processes-1)*n_samples_loc;
    }
    
    return n_samples_loc;
}

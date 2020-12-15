#include "definitions.hpp"

std::random_device rd;
std::mt19937_64 rng(rd());
std::uniform_int_distribution<int> binary(0,1);
std::uniform_real_distribution<double> unif(0.0,1.0);

vec4D flatten(mat2D M){
    
    vec4D v;
    v.head(2) = M.col(0);
    v.tail(2) = M.col(1);
    return v;
}

mat2D unflatten(vec4D v){
    
    mat2D M;
    M.col(0) = v.head(2);
    M.col(1) = v.tail(2);
    return M;
}


void display_spin_up(){
    
    printf("\033[34m%2s\033[0m", "O");
}

void display_spin_down(){
    
    printf("\033[37m%2s\033[0m", "X");
}

void print_bold(std::string message){
    
    printf("\033[1;37m%2s\033[0m", message.c_str());
}

void print_error(std::string message){
    
    printf("\033[1;31mERROR: %2s\033[0m\n", message.c_str());
}

void print_vec2D(vec2D v){
    
    printf("(%.7lf, %.7lf)\n", v(0), v(1));
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

std::string get_rounded_str(double num){
    
    std::stringstream ss;
    ss << std::setprecision(5) << num;
    std::string str = ss.str();
    
    return str;
}

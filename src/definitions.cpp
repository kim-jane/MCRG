#include "definitions.hpp"

std::random_device rd;
std::mt19937_64 rng(rd());
std::uniform_int_distribution<int> binary(0,1);
std::uniform_real_distribution<double> unif(0.0,1.0);


void display_spin_up(){
    
    printf("\033[37m%2s\033[0m", "O");
}

void display_spin_down(){
    
    printf("\033[34m%2s\033[0m", "O");
}

void print_bold(std::string message){
    
    printf("\033[1;37m%2s\033[0m", message.c_str());
}

bool print(int i){
    
    bool print = false;
    
    if(i < 10){
        print = true;
    }
    else if(i < 100 && i%10 == 0){
        print = true;
    }
    else if(i < 1000 && i%100 == 0){
        print = true;
    }
    else if(i < 10000 && i%1000 == 0){
        print = true;
    }
    else if(i < 100000 && i%10000 == 0){
        print = true;
    }
    else if(i%100000 == 0){
        print = true;
    }
    
    return print;
}

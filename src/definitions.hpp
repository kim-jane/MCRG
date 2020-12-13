#pragma once
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include "Eigen/Dense"

using namespace Eigen;


// linear algebra definitions
using mat = Matrix<double, Dynamic, Dynamic, ColMajor>;
using imat = Matrix<int, Dynamic, Dynamic, ColMajor>;
using mat2D = Matrix<double, 2, 2>;

using vec = Matrix<double, Dynamic, 1>;
using ivec = Matrix<int, Dynamic, 1>;
using vec2D = Matrix<double, 2, 1>;


// random number generators
extern std::random_device rd;
extern std::mt19937_64 rng;
extern std::uniform_int_distribution<int> binary;
extern std::uniform_real_distribution<double> unif;

inline auto rand_spin(){
    return 2*binary(rng)-1;
}

inline auto rand_unif(){
    return unif(rng);
}

extern void display_spin_up();
extern void display_spin_down();
extern void display_cluster_up();
extern void display_cluster_down();
extern void print_bold(std::string message);
extern void print_error(std::string message);
extern void print_vec2D(vec2D v);
extern bool write_iter(int i);
extern std::string get_rounded_str(double num);

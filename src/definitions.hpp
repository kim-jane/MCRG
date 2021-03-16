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

extern vec flatten(mat M);
extern mat unflatten(vec v);


// random number generators
extern std::random_device rd;
extern std::mt19937_64 rng;
extern std::uniform_int_distribution<int> binary;
extern std::uniform_real_distribution<double> unif;

extern std::normal_distribution<double> narrow_norm;

inline auto rand_spin(){
    return 2*binary(rng)-1;
}

inline auto rand_weight(){
    return narrow_norm(rng);
}

inline auto rand_unif(){
    return unif(rng);
}

extern void display_spin_up();
extern void display_spin_down();
extern bool write_iter(int i);
extern std::string get_rounded_str(double num, int precision);
extern int split_samples(int rank, int n_processes, int n_samples);

// element-wise square
extern vec square(vec v);

// element-wise squareroot
extern vec root(vec v);

// element-wise add double
extern vec add(vec v, double epsilon);
 
// element-wise divide
extern vec divide(vec u, vec v);

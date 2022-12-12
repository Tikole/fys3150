#pragma once

#include "slitex.hpp"

#define ARMA_USE_HDF5
#include <armadillo>

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>

using arma::cx_double;
using arma::cx_vec;
using arma::cx_mat;
using arma::sp_cx_mat;
using arma::vec;
using arma::ivec;
using arma::mat;

/* Make A and B matrices defining the time evolution of the wave packet*/
void make_A_and_B(const mat& V, double h, double Dt, sp_cx_mat& A_out, sp_cx_mat& B_out);
/* Matrix discretization of the potential */
void make_V(const SEParams& params, mat& V_out);
/* Initial state of the wave function */
void make_S0(const SEParams& params, cx_vec& S_out);
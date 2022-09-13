#pragma once

/* Supplies a function for solving tridiagonal matrix equations */

#include<armadillo>

void solve_tridiagonal(arma::vec& a, arma::vec& b, arma::vec& c, arma::vec& g);

void solve_special(arma::vec& g);

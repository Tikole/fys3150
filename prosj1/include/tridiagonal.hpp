#pragma once

/* Supplies a function for solving tridiagonal matrix equations */

#include<armadillo>

arma::vec solve_tridiagonal(arma::mat A, arma::vec g);

#pragma once

#include<armadillo>

#include<string>

const double pi = 3.14159265358979323846;

/* Makes 'output' into a tridiagonal matrix with subdiagonal elements all equal
to a, diagonal elements all equal to b, and superdiagonal elements all equal to
c. 'output' must be a NxN matrix with N > 1. A std::invalid_argument is thrown
if 'output' has wrong dimensions.*/
arma::mat& make_tridiag(double a, double b, double c, arma::mat& output);

/* Makes a dense symmetric matrix with dimensions of output, filled with random
numbers from a normal distribution. */
arma::mat& make_dense_symmetric(arma::mat& output);

/* Finds the largest magnitude element above the main diagonal. In a symmetric
matrix this is equivalent to finding the largest element off the main diagonal.
Returns the absolute value of this element, and writes its indexes i,j to
arguments i and j. */
double max_above_diagonal(const arma::mat& M, int& i, int& j);

/* Performs a similarity transform S^-1MS on the symmetric matrix M to eliminate element M(k,l),
where k!=l. Also performs the underlying rotation S to L.*/
void jacobi_rotate(arma::mat& M, arma::mat& L, int k, int l);

/*Solves for eigenvalues and eigenvectors of NxN matrix M. Eigenvalues outputed sorted ascendingly
in 'eigenvalues', with the eigenvector corresponding to the i'th eigenvalue, normalized, in the i'th
column of 'eigenvecs'. Returns number of similiarity transforms performed. Argument eps sets zero
threshold. */
int jacobi_method(arma::mat& M, arma::vec& eigenvalues, arma::mat& eigenvecs, double eps=1e-8);

void scaling_test(int N, std::string matrixtype);

/* Test functions. All return count of errors. Zeros means test passed.*/
int test_1(); /* Test function for make_tridiag() and arma::eig_sym().*/
int test_2(); // Tests max_above_diagonal()
int test_3(); // Test for jacobi_method() and jacobi_rotate()
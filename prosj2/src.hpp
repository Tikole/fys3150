#include<armadillo>

const double pi = 3.14159265358979323846;

/* Makes 'output' into a tridiagonal matrix with subdiagonal elements all equal
to a, diagonal elements all equal to b, and superdiagonal elements all equal to
c. 'output' must be a NxN matrix with N > 1. A std::invalid_argument is thrown
if 'output' has wrong dimensions.*/
arma::mat& make_tridiag(double a, double b, double c, arma::mat& output);

/* Finds the largest magnitude element above the main diagonal. In a symmetric
matrix this is equivalent to finding the largest element off the main diagonal.
Returns the value of this element, and writes its indexes row i,
column j to arguments i and j. */
double max_above_diagonal(const arma::mat& M, int& i, int& j);

/* Performs a similarity transform on the symmetric matrix M rotating away 
element M(k,l), where k!=l */
void jacobi_rotate(arma::mat& M, int k, int l);

/* Test functions. All return count of errors. Zeros means test passed.*/
int test_1(); /* Test function for make_tridiag() and arma::eig_sym().*/
int test_2(); // Tests max_above_diagonal()
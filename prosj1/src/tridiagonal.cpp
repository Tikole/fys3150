#include "tridiagonal.hpp"

#include<iostream>

arma::vec solve_tridiagonal(arma::mat A, arma::vec g) {
    /* Gaussian elimination assuming A is tridiagonal.
       All assignments of 0 or 1 to an element are not performed but implied,
       but implied by later operations. */
    const int n = A.n_rows;
    //The first row has no subdiagonal element, it only needs scaling.
    A(0,1) /= A(0,0);
    g(0) /= A(0,0;)
    // From top to bottom eliminate the subdiagonal element fromy every row
    // and scale to normalize pivot element.
    for (int i = 1; i < n; ++i) {
        g(i) = (g(i) - g(i-1)*A(i, i-1))/(A(i,i) - A(i-1,i)*A(i,i-1));
    }
    // Eliminate superdiagonal elements second bottommost to top.
    for (int i = n-2; i > 0; --i) {
        g(i) -= g(i-1)*A(i-1,i);
    }

    return g;
}
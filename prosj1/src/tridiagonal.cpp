#include "tridiagonal.hpp"

#include<iostream>

void solve_tridiagonal(arma::vec& a, arma::vec& b, arma::vec & c, arma::vec& g) {
    /* Gaussian elimination assuming A is a tridiagonal nxn matrix given by its
       subdiagonal, main diagonal, and superdiagonal elements. a, b, and c
       respectively. Solves the equation Av=g for v.
       b and g has length n. a and c must have length n-1.
       All assignments of 0 or 1 to an element are not performed,
       but implied by later operations.
       WIP: Might fail, the final result might not be a solution if no solution
       exists.*/
    // main diagonal element of row i is b(i), superdiagonal of row i is c(i),
    // and subdiagonal of row i is a(i-1)
    const int n = b.n_elem;
    //The first row has no subdiagonal element, it only needs scaling.
    c[0] /= b[0];
    g[0] /= b[0];
    // From top to bottom, except last row, eliminate the subdiagonal element fromy every row
    // and scale to normalize pivot element.
    for (int i = 1; i < n-1; ++i) {
        double scale_fac = b[i] - a[i-1]*c[i-1];
        g[i] = (g[i] - a[i-1]*g[i-1]) / scale_fac;
        c[i] /= scale_fac;
    }
    // Last row has no superdiagonal, must be handled specially to avoid out of bounds.
    int i = n - 1;
    g[i] = (g[i] - a[i-1]*g[i-1]) / (b[i] - a[i-1]*c[i-1]);
    // Eliminate superdiagonal elements second bottommost to top.
    for (int i = n-2; i >= 0; --i) {
        g[i] -= g[i+1]*c[i];
    }
}

void solve_special(arma::vec& g) {
    const int n = g.n_elem;
    const double a = -1;
    const double b = 2;
    const double c = -1;

    arma::vec cp(n-1);
    for (int i = 1; i < n; ++i) {
        cp[i-1] = -i/(i+1.0);
    }
    g[0] = g[0]/b;
    int i = 1;
    while(i < g.n_elem) {
        g[i] = (g[i]+g[i-1]) / (b - a*cp[i-1]);
        ++i;
    }
    i = n-2;
    while (i > -1) {
        g[i] -= cp[i]*g[i+1];
        --i;
    }
}
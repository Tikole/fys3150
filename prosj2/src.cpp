#include "src.hpp"

#include<cmath>
#include<iostream>
#include<string>

arma::mat& make_tridiag(double a, double b, double c, arma::mat& out) {
    if (out.n_rows != out.n_cols || out.n_rows < 2) {
        throw std::invalid_argument("A tridiagonal matrix must be square with"
                                    "N > 1");
    }
    out.zeros();
    out.diag(-1).fill(a);
    out.diag(0).fill(b);
    out.diag(1).fill(c);

    return out;
}

double max_above_diagonal(const arma::mat& M, int& I, int& J) {
    int n = M.n_cols;
    int m = M.n_rows;
    double max = 0.0;
    for (int j=1; j<n; ++j) {
        for (int i=0; (i < j) && (i < m); ++i) {
            if (std::abs(M(i,j)) > std::abs(max)) {
                max = M(i,j);
                I = i;
                J = j;
            }
        }
    }
    return max;
}

void jacobi_rotate(arma::mat& M, int k, int l) {
    const double tau = (M(l,l) - M(k,k))/(2*M(k,l));
    const double t = std::sqrt(1 + std::pow(tau,2)) - tau;
    const double c2 = 1/(1+std::pow(t,2));
    const double c = std::sqrt(c2);
    const double s = t*c;
    const double s2 = std::pow(s,2);
    const double N = M.n_rows;

    

    for (int i = 0; i < N; ++i) {
        Mik = M[i,k];
        Mil = M[i,l];
        Mkk = M[]
        if (i == k) {
            M[k,k] = Mik*c2 - 2*Mil*c*s + M[l,l]*s2;
            M[k,l] = 0.0;
        else if (i == l) {
            M[l,l] = M[l,l]*c2 + 2*M[k,l]*c*s + M[k,k]*s2;
        }
        }
    }
}

int test_1() {
    const int N = 6; // Size of matrix
    const double a = -1;
    const double d = 2;
    arma::mat M(N,N);
    make_tridiag(a, d, a, M);

    arma::vec eigval;
    arma::mat eigvec;
    bool b = arma::eig_sym(eigval, eigvec, M);

    eigval = arma::normalise(eigval);
    eigvec = arma::normalise(eigvec);

    arma::vec correct_eigval(N);
    arma::mat correct_eigvec(N,N);

    for (int i = 1; i <= N; ++i) {
        correct_eigval(i-1) = 2*(1 - std::cos((i*pi)/(N+1.0)));
        for (int j = 1; j <= N; ++j) {
            correct_eigvec(j-1,i-1) = std::sin((i*j*pi)/(N+1));
        }
    }

    correct_eigval = arma::normalise(correct_eigval);
    correct_eigvec = arma::normalise(correct_eigvec);

    int errors = 0;
    double tol = 1e-6;
    if (!arma::approx_equal(eigval, correct_eigval, "absdiff", tol)) {
        ++errors;
    }
    const arma::mat& A = eigvec;
    const arma::mat& B = correct_eigvec;
    for (int i = 0; i < A.n_cols; ++i) {
        if (!(arma::approx_equal(A.col(i), B.col(i), "absdiff", tol)) 
            && !(arma::approx_equal(A.col(i), -B.col(i), "absdiff", tol))) {
                ++errors;
            }
    }

    return errors;
}

int test_2() {
    arma::mat A{{1.0,  0.0,  0.0,  0.5},
                {0.0,  1.0, -0.7,  0.0},
                {0.0, -0.7,  1.0,  0.0},
                {0.5,  0.0,  0.0,  1.0}};
    
    double eps = 1e-6;
    int errors = 0;
    int i = 0;
    int j = 0;
    double d = max_above_diagonal(A, i, j);
    if (abs(d - A(1,2)) > eps)  {
        ++errors;
    }
    if (i != 1 || j != 2) {
        ++errors;
    }
    return errors;
}


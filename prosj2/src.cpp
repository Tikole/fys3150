#include "src.hpp"

#include<chrono>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<limits>
#include<string>

arma::mat& make_tridiag(double a, double b, double c, arma::mat& out) {
    if (out.n_rows != out.n_cols || out.n_rows < 2) {
        throw std::invalid_argument("A tridiagonal matrix must be square with"
                                    "N > 1");
    }
    out.zeros();
    out.diag(-1).fill(a); // Fill subdiagonal
    out.diag(0).fill(b); //  Fill main diagonal
    out.diag(1).fill(c); //  Fill superdiagonal

    return out;
}

arma::mat& make_dense_symmetric(arma::mat& out) {
    if (out.n_rows != out.n_cols || out.n_rows < 2) {
        throw std::invalid_argument("Matrix must be square, N>=2");
    }
    return out = arma::symmatu(arma::randu(out.n_rows, out.n_cols));
}

double max_above_diagonal(const arma::mat& M, int& I, int& J) {
    int n = M.n_cols;
    int m = M.n_rows;
    double max = 0.0;
    /* We find the largest element going column-wise, comparing all elements above the diagonal.*/
    for (int j=1; j<n; ++j) {
        for (int i=0; i<j && i<m; ++i) {
            if (std::abs(M.at(i,j)) > max) {
                max = std::abs(M.at(i,j));
                I = i;
                J = j;
            }
        }
    }
    return max;
}

void jacobi_rotate(arma::mat& M, arma::mat& L, int k, int l) {
    /* Applies the transformation S^TMS to M and SL to L. S is a transformation that rotates a
    column vector an amount theta in one plane such that element k,l of M is eliminated by the
    transformation S^TMS of M. */
    if (M.at(k,l)) { // If M_(k,l) is zero this is a null op.
                     //We must catch it to avoid division by zero.
        const double tau = (M.at(l,l) - M.at(k,k))/(2*M.at(k,l)); // Parameter that essentially
                                                                  // decides the amount of rotation.
        const double t = (tau > 0) // t is tan(theta)
                ? (-tau + std::sqrt(tau*tau + 1)) : (-tau - std::sqrt(tau*tau + 1));
        const double c2 = 1/(1+t*t); // cos^2(theta)
        const double c = std::sqrt(c2);
        const double s = t*c; // sin(theta)
        const double s2 = s*s;

        /* The transformation of these elements can be taken out of the loop over the rows. */
        const double Mkk = M.at(k,k);
        const double Mll = M.at(l,l);
        M.at(k,k) = c2*Mkk - 2*c*s*M.at(k,l) + s2*Mll;
        M.at(l,l) = c2*Mll + 2*c*s*M.at(k,l) + s2*Mkk;
        /* The rotation is determined precisely so that these will be zero, so no need to calculate. */
        M.at(k,l) = 0.0;
        M.at(l,k) = 0.0;

        const double N = M.n_rows;
        /* Apply transformation to all elements of column k, and l, except those handled above. Other
        columns are not touched by the transformation. */
        for (int i = 0; i < N; ++i) {
            if (i != k && i != l) {
                double Mik = M.at(i,k);
                double Mil = M.at(i,l);
                M.at(i,k) = c*Mik - s*Mil;
                M.at(k,i) = M.at(i,k);
                M.at(i,l) = c*Mil + s*Mik;
                M.at(l,i) = M.at(i,l);
            }
            /* Apply rotation column k and l of L.*/
            double Lik = L.at(i,k);
            double Lil = L.at(i,l);
            L.at(i,k) = c*Lik - s*Lil;
            L.at(i,l) = c*Lil + s*Lik;
        }
    }
}

void pm(const arma::mat& M) {
    M.print();
} 
void pv(const arma::vec& v) {
    v.print();
}

int jacobi_method(arma::mat& M, arma::vec& eigval, arma::mat& eigvec, double eps) {
    int m = M.n_rows;
    int n = M.n_cols;
    if (m != n || m < 2) {
        throw std::invalid_argument("Matrix must be square NxN, with N>=2");
    }
    eigvec = arma::eye(m, n);
    /* k and l will hold indexes of the element returned by max_above_diagonal(). d will
    hold this elements magnitude. */
    int k = -1;
    int l = -1;
    double d = -1.0;
    /* We'll keep track of the number of transformations with it, and set a limit to how
    many iterations we allow to prevent hanging. */
    int it = 0;
    int max_it = 1e6;
    /* Keep transforming until all non-diagonals are smaller than eps. */
    while(d = max_above_diagonal(M,k,l) > eps && it++ < max_it) {
        jacobi_rotate(M, eigvec, k, l);
    }
    eigvec = arma::normalise(eigvec);
    eigval.set_size(n);
    eigval = M.diag();
    /* Sorter egenverdier/-vektorer*/
    int ne = eigval.n_elem;
    for (int i = 0; i<ne-1; ++i) { // For hver i: finn det elementet blant i,i+1,... som er minst.
        double smallest = eigval(i);
        int smallest_index = i;
        for (int j = i+1; j<ne; ++j) {
            if (eigval.at(j) < smallest) {
                smallest = eigval.at(j);
                smallest_index = j;
            }
        }
        if (i != smallest_index){
            eigval.swap_rows(i,smallest_index);
            eigvec.swap_cols(i,smallest_index);
        }
    }

    return (it <= max_it) ? (it) : (-it); // If surpassed maximum allowed iterations return -it as sign of failure
}

void scaling_test(int N, std::string mt){
    if (mt != std::string("scaling") && mt != std::string("densescaling"))
        throw std::invalid_argument("Invalid matrix type");
    std::string type = (mt == std::string("scaling")) 
        ? (std::string("tridiagonal symmetric")) : (std::string("dense symmetric"));
    std::cout << "Using Jacobi method to calculate eigenvalues/-vectors\n"
              << "for NxN matrix of type: " << type << "\n---\n";
    arma::ivec Nl = arma::regspace<arma::ivec>(3, N); 
    arma::vec tl = arma::zeros(Nl.n_elem);
    arma::ivec sl = arma::zeros<arma::ivec>(Nl.n_elem);
    for (int i=0; i<Nl.n_elem; ++i) { 
        int n = Nl(i);
        std::string type = (mt == std::string("scaling"))
            ? ("tridiagonal symmetric") : ("dense symmetric");
        std::cout << "N = " << n << "\n";
        std::cout << "Solving...\n";
        arma::mat A(n,n);
        arma::vec eigval;
        arma::mat eigvec;
        if (mt == std::string("scaling"))
            make_tridiag(-1,2,-1,A);
        else
            make_dense_symmetric(A);

        auto t1 = std::chrono::high_resolution_clock::now();
        int it = jacobi_method(A,eigval,eigvec);
        auto t2 = std::chrono::high_resolution_clock::now();
        double Dt = std::chrono::duration<double>(t2 - t1).count();

        if (it > 0) {
            std::cout << "Solved in " << Dt << " seconds, using " << it << " transformations.\n\n";
            tl(i) = Dt;
            sl(i) = it;
        }
        else {
            std::cout << "Timed out after " << Dt << " seconds, hitting limit of "
                      << -it-1 << " transformations\n\n";
            break; // We'll assume the rest will fail as well.
        }
    }
    std::fstream file;
    std::string filename =
        (mt == std::string("scaling")) ? ("scaling_tridiag.txt") : ("scaling_dense.txt");
    file.open(filename, std::ios::out | std::ios::trunc);
    int d = std::numeric_limits<double>::digits10; // Use all the precision we have
    file << std::scientific << std::setprecision(d);
    for (int i=0; i<Nl.n_elem; ++i){
        if (sl(i) > 0)
            file << Nl(i) << ' ' << sl(i) << ' ' << tl(i) << std::endl;
    }
    file.close();

}

int test_1() {
    // Make a test matrix.
    const int N = 6;
    const double a = -1;
    const double d = 2;
    arma::mat M(N,N);
    make_tridiag(a, d, a, M);

    // Calculate eigenvalues/-vectors with armadillo
    arma::vec eigval;
    arma::mat eigvec;
    bool b = arma::eig_sym(eigval, eigvec, M);

    eigvec = arma::normalise(eigvec);

    /* Calculate the correct results using an anlytical formula.*/
    arma::vec correct_eigval(N);
    arma::mat correct_eigvec(N,N);

    for (int i = 1; i <= N; ++i) {
        correct_eigval(i-1) = 2*(1 - std::cos((i*pi)/(N+1.0)));
        for (int j = 1; j <= N; ++j) {
            correct_eigvec(j-1,i-1) = std::sin((i*j*pi)/(N+1));
        }
    }
    correct_eigvec = arma::normalise(correct_eigvec);

    /* Compare eigenvalues. */
    int errors = 0;
    double tol = 1e-6;
    if (!arma::approx_equal(eigval, correct_eigval, "absdiff", tol)) {
        ++errors;
    }
    /* Compare eigenvectors. */
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
    /* Test matrix. */
    arma::mat A{{1.0,  0.0,  0.0,  0.5},
                {0.0,  1.0, -0.7,  0.0},
                {0.0, -0.7,  1.0,  0.0},
                {0.5,  0.0,  0.0,  1.0}};
    
    double eps = 1e-8;
    int errors = 0;
    int i = 0;
    int j = 0;
    /* The correct result is obviously i=1, j=2, d=0.7.*/
    double d = max_above_diagonal(A, i, j);
    if (d - std::abs(A(i,j)) > eps)  {
        ++errors;
    }
    return errors;
}

int test_3() {
    /*Test matrix*/
    const int N = 6;
    arma::mat M(N,N);
    make_tridiag(-1,2,-1,M);

    /* Calculate using our implementation of the jacobi method. */
    arma::mat eigvecs;
    arma::vec eigvals;
    arma::mat M2(M);
    jacobi_method(M, eigvals, eigvecs);

    /* Correct values using analytical formula. */
    arma::vec correct_eigvals(N);
    arma::mat correct_eigvecs(N,N);
    for (int i = 1; i <= N; ++i) {
        correct_eigvals(i-1) = 2*(1 - std::cos((i*pi)/(N+1.0)));
        for (int j = 1; j <= N; ++j) {
            correct_eigvecs(j-1,i-1) = std::sin((i*j*pi)/(N+1));
        }
    }
    correct_eigvecs = arma::normalise(correct_eigvecs);

    // arma::vec correct_eigvals_2;
    // arma::mat correct_eigvecs_2;
    // arma::eig_sym(correct_eigvals_2, correct_eigvecs_2, M2);
    // correct_eigvecs_2 = arma::normalise(correct_eigvecs_2);

    // std::cout << "Calculated eigenvalues:\n";
    // eigvals.print();
    // std::cout << "\n---\n";
    // std::cout << "Correct eigenvalues:\n";
    // correct_eigvals.print();
    // std::cout << "\n---\n";
    // std::cout << "Arma eigenvalues:\n";
    // correct_eigvals_2.print();
    // std::cout << "\n---\n";
    // std::cout << "Calculated eigenvectors:\n";
    // eigvecs.print();
    // std::cout << "\n---\n";
    // std::cout << "Correct eigenvectors:\n";
    // correct_eigvecs.print();
    // std::cout << "\n---\n";
    // std::cout << "Arma eigenvectors:\n";
    // correct_eigvecs_2.print();
    // std::cout << "\n---\n";

    /* Compare eigenvalues. */
    int errors = 0;
    double tol = 1e-6;
    if (!arma::approx_equal(eigvals, correct_eigvals, "absdiff", tol)) {
        ++errors;
    }
    /* Compare eigenvectors. */
    const arma::mat& A = eigvecs;
    const arma::mat& B = correct_eigvecs;
    for (int i = 0; i < A.n_cols; ++i) {
        if (!(arma::approx_equal(A.col(i), B.col(i), "absdiff", tol)) 
            && !(arma::approx_equal(A.col(i), -B.col(i), "absdiff", tol))) {
                ++errors;
            }
    }

    return errors;
}


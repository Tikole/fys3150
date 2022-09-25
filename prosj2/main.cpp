#include "src.hpp"

#include<armadillo>

#include<fstream>
#include<iomanip>
#include<iostream>
#include<limits>
#include<string>
#include<vector>


void print_usage() {
    std::cout << "NOT YET IMPLEMENTED" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc == 1) { // Only one arg: print usage and quit
        print_usage();
        return 0;
    }
    else{
        std::string command = std::string(argv[1]);
        if (command == std::string("test")) { // RUn all tests, print error counts
            std::cout << "Running tests\n---\n"
                      << "Test 1...\n"
                      << test_1() << " errors.\n";
            std::cout << "Test 2...\n"
                      << test_2() << " errors.\n";
            std::cout << "Test 3...\n"
                      << test_3() << " errors.\n";

        }
        else if (command == std::string("scaling")) { 
            /* Test scaling of our implementation of the jacobi method.
            After the 'scaling' command there should be a single integer
            argument n > 2. Jacobi method is used on matrixes of dimensions
            3x3, 4,4,..., NxN.*/
            if (argc != 3) {
                std::cout << "The 'scaling' command should be followed by a single integer argument.";
            }
            else {
                int n;
                try {
                    n = std::stoi(argv[2]);
                }
                catch(...) {
                    std::cout << "The 'scaling' command should be followed by a single integer argument.";
                    return 0;
                }
                scaling_test(n, TRIDIAGONAL);
                scaling_test(n, DENSE);
            }
        }
        else if (command == std::string("solve")) {
            if (argc != 3) {
                std::cout << "The 'solve' command should be followed by a single integer argument";
            }
            else {
                int steps;
                try {
                    steps = std::stoi(argv[2]);
                }
                catch(...) {
                    std::cout << "The 'solve' command should be followed by a single integer argument";
                    return 0;
                }
                int N = steps - 1;
                double h2 = 1.0/(steps*steps);
                arma::mat A(N, N);
                make_tridiag(-1/h2, 2/h2, -1/h2, A);
                arma::vec eigenval;
                arma::mat eigenvec;
                int it = jacobi_method(A, eigenval, eigenvec);
                if (it < 0) {
                    std::cout << "Solve failed after hitting limit of " << -it << " transformations.";
                    return 0;
                }
                eigenvec.insert_rows(0, arma::rowvec(eigenvec.n_cols, arma::fill::zeros));
                eigenvec.insert_rows(eigenvec.n_rows, arma::rowvec(eigenvec.n_cols, arma::fill::zeros));

                arma::vec analytical_eigenval;
                arma::mat analytical_eigenvec;
                analytical_solution(N, analytical_eigenval, analytical_eigenvec);
                analytical_eigenvec.insert_rows(0,
                                                arma::rowvec(analytical_eigenvec.n_cols, arma::fill::zeros)
                );
                analytical_eigenvec.insert_rows(analytical_eigenvec.n_rows,
                                                arma::rowvec(analytical_eigenvec.n_cols, arma::fill::zeros)
                );
                std::string filename = "solution_" + std::to_string(steps) + ".txt";
                std::fstream file(filename, std::ios::out | std::ios::trunc);
                int d = std::numeric_limits<double>::digits10; // Use all possible precision.
                file << std::scientific << std::setprecision(d);
                for (int i=0; i<eigenvec.n_rows; ++i) {
                    double x = i*std::sqrt(h2);
                    file << x << ' ' << eigenvec(i,0) << ' ' << eigenvec(i,1) << ' '
                         << eigenvec(i,2) << ' ' << analytical_eigenvec(i,0) << ' '
                         << ' ' << analytical_eigenvec(i,1) << ' '
                         << analytical_eigenvec(i,2) << std::endl;
                }
            }
        }
        else {
            std::cout << "Unrecognized command: '" << std::string(argv[1]) << "'. For help provide no arguments.\n";
        }
    }

    return 0;
}
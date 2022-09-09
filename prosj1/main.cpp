#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<limits>

#include "tridiagonal.hpp"

double u(double x){
    return 1.0 - (1 - std::exp(-10.0))*x - std::exp(-10.0*x);
}

double f(double x){
    return 100*std::exp(-10.0*x)
}

void print_usage() {
    std::cout << "Evaluates u(x) defined by the poisson equation:\n"
    "    -d^2u/dx^2 = f(x), where f(x)=100e^(-10x), u's domain is [0,1], "
    "and u(0) = u(1) = 0\n"
    "at the number of evaluation points provided as arguments, using both the"
    "analytical solution as well as solving the discretized equation using the"
    "\"tridiagonal.hpp\" tridiagonal matrix equation solver."
    "Output is provided in .txt files with two whitespace separated columns of x- and y-values\n"
    "A single \"analytical_output.txt\" is outputted, and a number of"
    "\"discretized_output_<eval_points>.txt\".\n"
    "A sequence of integer arguments specifying number of evaluations points must be provided.\n"
    "The first argument is used for the analytical solution, the rest are evaluated discretized.";
}

int main(int argc, char* argv[]){
    /* Check arguments */
    int n_analytical = 0;
    std::vector<int> n_discreet;
    try {
        if (argc == 0) {
            print_usage();
        }
        else if (argc == 1) {
            n_analytical = std::stoi(argv[1]);
        }
        else {
            for (int i=2; i<argc; ++i) {
                n_discreet.emplace_back(std::stoi(argv[i]));
            }
        }
    }
    catch (...) {
        std::cout << "Invalid argument. Arguments should be a list of integers.";
        return 0;
    }

    /* Evaluate function between, and including, x0 and x1. */
    const double x0 = 0.0;
    const double x1 = 1.0;
    double dx = (x1 - x0)/(n_analytical-1);
    std::vector<double> X(n_analytical);
    int i = 0;
    while (i < X.size()-1) {
        X[i] = x0 + i*dx;
        ++i;
    }
    X[i] = x1;

    std::vector<double> Y(n);
    for (auto i = 0; i < X.size(); ++i) {
        Y[i] = u(X[i]);
    }

    /* Write numbers with all the precision allowed by the double type in scientific format. */
    std::fstream output_file("analytical_output.txt", std::ios::trunc | std::ios::out);
    output_file << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
    for (int i = 0; i < X.size(); ++i){
        output_file << X[i] << " " << Y[i] << std::endl;
    }
    output_file.close()
    
    /* Solve discretized equation*/
    for (auto it = n_discreet.begin(); it != n_discreet.end(); ++it) {
        n = *it; // Number of elements of solution
        /* Set up coefficient matrix. */
        arma::mat A(n-2, n-2);
        A(0,0) = 2;
        A(0,1) = -1;
        for (int i = 0; i < A.n_rows; ++i) {
            A(i,i) = 2;
            A(i,i-1) = A(i,i+1) = -1;
        }
        A(A.n_rows-1, A.n_cols-2) = -1;
        A(A.n_rows-1, A.n_cols-1) = 2;

        /* Set up g. Note that since the boundary conditions are zero we dont
           have to treat the first and last element specially. */
        arma::vec g(n-2);
        arma::vec X(n-2)
        double dx = (x1 - x0)/(n-1)
        for (int i = 0; i < g.n_rows; ++i) {
            X(i) = x0 + i*dx; 
            g(i) = -f(X(i))/std::pow(dx,2); 
        }

        arma::vec solution = solve_tridiagonal(A, g)
        std::fstream output(
            "discretized_output_" + to_string(n) + ".txt",
            std::ios::trunc | std::ios::out
        );
        output << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
        for (int i = 0; i < X.size(); ++i) {
            output << X[i] << " " << solution[i] << std::endl;
        }
    }

    return 0;
}
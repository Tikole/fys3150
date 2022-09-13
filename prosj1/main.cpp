#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<limits>
#include<stdexcept>

#include "tridiagonal.hpp"

void print_usage() {
    std::cout << "Evaluates u(x) defined by the poisson equation:\n"
    "    -d^2u/dx^2 = f(x), where f(x)=100e^(-10x), u's domain is [0,1], "
    "and u(0) = u(1) = 0\n"
    "at the number of evaluation points provided as arguments, using both "
    "the analytical solution\n"
    "as well as solving the discretized equation using the \"tridiagonal.hpp\""
    "special and general\n"
    "tridiagonal matrix equation solver. Output is provided in .txt files "
    "with four whitespace\n"
    "separated columns of x-value, analytical y-values, numerical y_values "
    "using a general solver,\n"
    "and numerical y-values using a special solver. Files are named "
    "\"output_<eval_points>.txt.\".\n"
    "A sequence of integer arguments specifying number of evaluations points "
    "must be provided.\n";
}

// Limits of domain.
const double x0 = 0;
const double x1 = 1;

double u(double x){
    return 1.0 - (1 - std::exp(-10.0))*x - std::exp(-10.0*x);
}

double f(double x){
    return 100*std::exp(-10.0*x);
}

void evaluate_u(arma::vec& output_x, arma::vec& output_y) {
    /* Evaluate u(x) at the number of points given by the length of output_x
    and output_y, who must have the same length. Initial contents of output_x
    and output_y are discarded. */
    int n = output_x.n_elem;
    if (n != output_y.n_elem) {
        throw std::invalid_argument(
            "output_x and output_y must have same length."
        );
    }
    double dx = (x1 - x0)/(n-1);
    for (int i = 0; i < n; ++i) {
        double x = x0 + i*dx;
        output_x(i) = x;
        output_y(i) = u(x);
    }
}

void make_tridiagonal_signature(arma::mat& A, double a, double b, double c) {
    /* Makes A into a tridiagonal matrix with a in all the subdiagonal, c in 
    the superdiagonal elements, b in all the main diagonal elements.*/
    int n = A.n_rows;
    if (n != A.n_cols || n < 2) {
        throw std::invalid_argument (
            "Can't make a  tridiagonal matrix with these dimensions."
        );
    }
    A(0,0) = b;
    A(0,1) = c;
    for (int i = 1; i < n-1; ++i) {
        A(i,i-1) = a;
        A(i,i) = b;
        A(i,i+1) = c;
    }
    A(n-1, n-2) = a;
    A(n-1, n-1) = b;
}

void write_file(int n,
                const arma::vec& X,
                const arma::vec& Ya,
                const arma::vec& Yg,
                const arma::vec& Ys)
{   
    /* Writes evaluation of u into a file named output_<n>.txt*/
    std::fstream output(
        "output_" + std::to_string(n) + ".txt",
        std::ios::trunc | std::ios::out
    );
    int d = std::numeric_limits<double>::digits10; // Use all the precision we have
    output << std::scientific << std::setprecision(d);
    int i = 0;
    output << X(i) << ' ' << Ya(i) << ' ' << Ya(i) << ' ' << Ya(i) << std::endl;
    while (++i < Yg.n_elem) {
        output << X(i) << ' ' << Ya(i) << ' ' << Yg(i-1) << ' ' << Ys(i-1)
               << std::endl;
    }
    output << X(i) << ' ' << Ya(i) << ' ' << Ya(i) << ' ' << Ya(i) << std::endl;
}

int main(int argc, char* argv[]){
    /* Check arguments */
    if (argc == 1) {
        print_usage();
        return 0;
    }

    /* Create a list of numbers of evaluation points*/
    std::vector<int> ln_eval_points(argc - 1);
    for (int i = 1; i < argc; ++i) {
        ln_eval_points.at(i-1) = std::stoi(argv[i]);
    }

    /* Evaluate, solve and write*/
    for (auto it = ln_eval_points.begin(); it != ln_eval_points.end(); ++it) {
        int n = *it; // Length of solution
        arma::vec X(n); // X-values to evaluate

        /*Analytical solution*/
        arma::vec Ya(n);
        evaluate_u(X, Ya);
        double dx = X(1) - X(0);

        arma::vec g(n-2);
        for (int i = 0; i < g.n_rows; ++i) {
            g(i) = f(X(i+1))*std::pow(dx, 2);
        }

        /* Solve using general tridiagonal solver*/
        arma::vec a(n-3);
        a.fill(-1.0);
        arma::vec b(n-2);
        b.fill(2.0);
        arma::vec c(n-3);
        c.fill(-1.0);
        arma::vec Yg(g);
        solve_tridiagonal(a, b, c, Yg);

        /* Solve using special solver*/
        solve_special(g);

        /*Write results to file*/
        write_file(n, X, Ya, Yg, g);
    }

    return 0;
}
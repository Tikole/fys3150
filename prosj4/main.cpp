#include "src.hpp"

#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <string>

void print_usage() {
    std::iostream << "NOT YET IMPLEMENTED\n";
}

int main(int argc, char* argv[]) {
    /* Handle arguments, set up parameters structure, call solver. */
    Parameters P;
    /* Handle positional arguments, all of which should be positive integers
    giving lattice sizes. */
    int i = 0;
    while (++i < argc) {
        // Found named argument, break
        if (argv[i][0] == '-')
            break;
        else 
            P.Lvec.insert_rows(i-1, std::stoi(argv[i]));
    }
    /* Named arguments */
    while (i < argc) {
        std::string s(argv[i]); 
        // Temperature as float, float, int specifying a 'linspace'
        if (s == "-Tspace") {
            double T0 = std::stod(argv[++i]);
            double T1 = std::stod(argv[++i]);
            int steps = std::stoi(argv[++i]);
            P.Tvec = arma::linspace(T0, T1, steps);
            ++i;
        }
        // Temperature as a single float
        else if (s == "-T") {
            double T = std::stod(argv[++i]);
            P.Tvec = arma::vec{T};
            ++i;
        }
        else if (s == "-random") {
            P.randomize = true;
        }   
    }

    multirun(P);

    return 0;
}

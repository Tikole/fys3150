#include "src.hpp"

#include <iostream>
#include <sstream>
#include <string>

/* Solves the one particle analytical case with both FE and RK4.
Solve for time t with n steps. */
void solve_analytical_case(double t, int n);
void solve_twoparticle_case(double t, int n, bool b);
void print_usage();

int main(int argc, char* argv[]) {
    if (argc == 1) {
        print_usage();
        return 0;
    }
    
    std::string command = std::string(argv[1]);
    if (command == "analytic") {
        /* Generate analytical case data */
        double t;
        int n;
        try {
            t = std::stod(argv[2]);
            n = std::stoi(argv[3]);
        }
        catch (...) {
            std::cout << "'analytic' command must be followed by 2 whitespace separated arguments:\n"
                      << "time and number of steps: ex. 'analytic 50.0 4000'\n";
            return 0;

        }
        solve_analytical_case(t, n);
    }
    else if (command == "twoparticle") {
        /* Generate data on case of two particles */
        double t;
        int n;
        int b;
        try {
            t = std::stod(argv[2]);
            n = std::stoi(argv[3]);
            b = std::stoi(argv[4]);
        }
        catch (...) {
            std::cout << "'twoparticle' command must be followed by 3 whitespace separated arguments:\n"
                      << "time, number of steps and coulomb interactions flag ex. 'twoparticle 50 4000 1'";
        }
        solve_twoparticle_case(t, n, b);
    }
    else {
        std::cout << "Invalid command. Run argumentless for help.";
    }
    return 0;
}

void print_usage() {
    std::cout << "NOT YET IMPLEMENTED" << std::endl;
}

void solve_analytical_case(double t, int n) {
    PenningTrap pt;
    PT_states FE = pt.advance_forward_euler(t, n);
    PT_states RK4 = pt.advance_runge_kutta_4(t, n);
    std::stringstream s1;
    s1 << "Analytical case t=" << (int)t << " us solved with forward Euler using "
       << n << " steps.";
    write_to_file(FE, "FE", "ana", s1.str());
    std::stringstream s2;
    s2 << "Analytical case t=" << (int)t << " us solved with Runge Kutta 4 using "
       << n << " steps.";
    write_to_file(RK4, "RK4", "ana", s2.str());
}

void solve_twoparticle_case(double t, int n, bool b) {
    Particles pv({
        Particle(Ca_mass, Ca_charge, {20.0, 0.0, 20.0}, {0.0, 25.0, 0.0}),
        Particle(Ca_mass, Ca_charge, {25.0, 25.0, 0.0}, {0.0, 40.0, 5.0})
    });
    PenningTrap pt(def_B, def_V, def_d, pv, b);
    PT_states twop = pt.advance_runge_kutta_4(t,n);
    std::stringstream s;
    s << "Two particles, " << ((b) ? ("interacting") : ("non-interacting"))
        << ", " << "t=" << t << " us solved with Runge Kutta 4 using " << n
        << " steps.";
    write_to_file(twop, "RK4", "twop", s.str());
}
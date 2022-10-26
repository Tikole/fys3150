#include "src.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

/* ***Solves cases for time t[us] and number of steps n:*** */
// One particle analytical case with both FE and RK4.
void solve_analytical_case(double t, int n);
// Two particles with/without coulomb interactions with RK4.
void solve_twoparticle_case(double t, int n, bool b);
// Of N particles with random initial velocity/position sampled from a normal
// distribution, in an electric field fluctuating sinusoidally, count how many
// are left.
void solve_nrand(double t, int n, bool b, int N, arma::vec A, arma::vec f);

void print_usage();

int main(int argc, char* argv[]) {
    if (argc == 1) {
        print_usage();
        return 0;
    }
    
    std::string command = std::string(argv[1]);
    if (command == "analytic") {
        /* Generate analytical case data */
        double t; // time
        int n; // steps
        try {
            t = std::stod(argv[2]);
            n = std::stoi(argv[3]);
        }
        catch (...) {
            std::cout << "Invalid arguments for 'analytic' command must be followed by 2 whitespace separated arguments:\n"
                      << "time and number of steps: ex. 'analytic 50.0 4000'\n";
            return 0;

        }
        solve_analytical_case(t, n);
    }
    else if (command == "twoparticle") {
        /* Generate data on case of two particles */
        double t; // time
        int n; // steps
        int b; // coulomb interactions?
        try {
            t = std::stod(argv[2]);
            n = std::stoi(argv[3]);
            b = std::stoi(argv[4]);
        }
        catch (...) {
            std::cout << "'twoparticle' command must be followed by 3 whitespace separated arguments:\n"
                      << "time, number of steps and coulomb interactions flag ex. 'twoparticle 50 4000 1'\n";
        }
        solve_twoparticle_case(t, n, b);
    }
    else if (command == "nrand") {
        double t; // time
        int n; // time steps
        int b; // coulomb interactions?
        int N; // Number of particles
        arma::vec f_vec; // frequencies
        arma::vec A_vec; // amplitudes
        try {
            t = std::stod(argv[2]);
            n = std::stoi(argv[3]);
            b = std::stoi(argv[4]);
            N = std::stoi(argv[5]);
            double f0 = std::stod(argv[6]);
            double f1 = std::stod(argv[7]);
            int nf = std::stod(argv[8]);
            f_vec = arma::linspace(f0, f1, nf);
            int nA = argc - 9;
            if (nA < 1)
                throw 0;
            A_vec.zeros(nA);
            for (int i = 9; i < argc; ++i) {
                A_vec(i-9) = std::stod(argv[i]);
            }
            A_vec = arma::sort(A_vec);
        }
        catch(...) {
            std::cout << "Invalid arguments for 'nrand'\n";
        }
        solve_nrand(t, n, b, N, A_vec, f_vec);
    }
    else {
        std::cout << "Invalid command. Run argumentless for help.\n";
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

void solve_nrand(double t, int n, bool b, int N, arma::vec A,
    arma::vec f)
{
    arma::imat M(A.n_elem, f.n_elem);
    for (int i=0; i<A.n_elem; ++i)
        for (int j=0; j<f.n_elem; ++j) {
            std::cout << i << '/' << A.n_elem << ' ' << j << '/' << f.n_elem << std::endl;
            // Make 100 new random particles.
            Particles pv(N);
            for (int i = 0; i<N; ++i)
                pv[i] = Particle{Ca_mass, Ca_charge, 0.1*def_d, 0.1*def_d};
            PenningTrap pt(def_B, def_V, f[j], A[i], 0.0, def_d, pv, b);
            PenningTrap end_state = pt.advance_runge_kutta_4(t, n).back();
            M(i,j) = end_state.n_trapped();
        }
    /* Write to file*/
    std::string filename = "nrand_" + std::to_string(int(b)) + ".txt";
    std::fstream file(filename, std::ios::trunc | std::ios::out);
    for (int i = 0; i < A.n_elem; ++i)
        file << A(i) << ' ';
    file << std::endl;
    for (int i=0; i<f.n_elem; ++i)
        file << f(i) << ' ';
    file << std::endl;
    for (int i=0; i<M.n_cols; ++i) {
        for (int j=0; j<M.n_rows; ++j) {
            file << M(j,i) << ' ';
        }
        file << std::endl;
    }
}
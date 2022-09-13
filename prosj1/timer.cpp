#include<chrono>
#include<vector>
#include<iostream>

#include<armadillo>

#include "tridiagonal.hpp"

int main() {
    /* n_steps to test at. */
    const std::vector<int> nsteps{100, 1000, 10000, 100000, 1000000};

    /*Number of repetitions at each n_steps*/
    const int nr = 10;

    /* Holds the times measured*/
    arma::vec tg(nr);
    arma::vec ts(nr);

    for (int i = 0; i < nsteps.size(); ++i) {
        int n = nsteps[i];
        std::cout << "For n = " << n << ':' << std::endl;
        for (int j = 0; j < nr; ++j) {
            // Both methods should be solving the same equation.
            arma::vec g1(n);
            g1.fill(0.01234);
            arma::vec g2(g1);

            // initialize the tridiagonal matrix
            arma::vec a(n-1);
            a.fill(-1.0);
            arma::vec b(n);
            b.fill(2.0);
            arma::vec c(n-1);

            // Measure the time used by the general solver.
            auto t1_g = std::chrono::high_resolution_clock::now();
            solve_tridiagonal(a,b,c,g1);
            auto t2_g = std::chrono::high_resolution_clock::now();
            double delta_tg = std::chrono::duration<double>(t2_g - t1_g).count();
            tg[j] = delta_tg;

            // Measure the used by the special solver.
            auto t1_s = std::chrono::high_resolution_clock::now();
            solve_special(g2);
            auto t2_s = std::chrono::high_resolution_clock::now();
            double delta_ts = std::chrono::duration<double>(t2_s - t1_s).count();
            ts[j] = delta_ts;
        }
        std::cout << "    General solver:\n"
                  << "        avg: " << arma::mean(tg) << "s\n"
                  << "        std.dev.: " << arma::stddev(tg) << "\n"
                  << "    Special solver:\n"
                  << "        avg: " << arma::mean(ts) << "s\n"
                  << "        std.dev.: " << arma::stddev(ts) << "\n---\n";
    }
}
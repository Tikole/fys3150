#include "src.hpp"

#include <armadillo>

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>

void print_usage() {
    std::cout << "NOT YET IMPLEMENTED" << std::endl;
}

int main(int argc, char* argv[]) {
    /* Handle arguments, set up parameters structure, call solver. */
    if (argc < 2) {
        std::cout << "No arguments provided. (-h for help)" << std::endl;
        return 0;
    }
    Parameters P;
    if (std::string(argv[1]) == "-defaults") {
        std::cout << "-random" << std::endl
                  << "-runs " << P.runs << std::endl
                  << "-max " << P.max_chunks << std::endl
                  << "-target " << P.target_chunks << std::endl
                  << "-chunk " << P.chunk_sz << std::endl
                  << "-sample " << P.sample_frequency << std::endl
                  << "-eps " << P.epsilon << std::endl;
        return 0;
    }
    if (std::string(argv[1]) == "-help" || std::string(argv[1]) == "-h") {
        print_usage();
        return 0;
    }
    try {
        /* Handle positional arguments, a single int providing L */
        int i = 1;
        P.L = std::stoi(argv[i++]);
        
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
                ++i;
            }
            else if (s == "-ordered") {
                P.randomize = false;
                ++i;
            }
            else if (s == "-runs") {
                P.runs = std::stoi(argv[++i]);
                ++i;
            }
            else if (s == "-max") {
                P.max_chunks = std::stoi(argv[++i]);
                ++i;
            }
            else if (s == "-target") {
                P.target_chunks = std::stoi(argv[++i]);
                ++i;
            }
            else if (s == "-chunk") {
                P.chunk_sz = std::stoi(argv[++i]);
                ++i;
            }
            else if (s == "-sample") {
                P.sample_frequency = std::stoi(argv[++i]);
                ++i;
            }
            else if (s == "-eps") {
                P.epsilon = std::stod(argv[++i]);
                ++i;
            }
            else
                ++i;

        }
    }
    catch (...) {
        std::cout << "ERROR: Ill-formed parameters" << std::endl;
    }

    try {
        auto t1_g = std::chrono::high_resolution_clock::now();
        multirun(P);
        auto t2_g = std::chrono::high_resolution_clock::now();
        double delta_tg = std::chrono::duration<double>(t2_g - t1_g).count();
        std::cout << "Time elapsed: " << delta_tg << 's' << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
    }

    return 0;
}

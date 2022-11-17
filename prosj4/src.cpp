#include "src.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

/* Save results to this directory */
const char* data_directory= "./data";

/* Data structure serving to pass arguments to a run and return results thereof*/
struct Results {
    int L;
    double T;
    int equilibrium_index;
    int max_cycles;
    int target_cycles;
    int target_chunks;
    int chunk_sz;
    int sample_frequency;
    double epsilon;
    arma::ivec E;
    arma::ivec M;
    arm::vec E_chunks;
    
    void initialize(const Parameters& P, int L_index, int T_index);
};

void Results::initialize(const Parameters& P, int L_index, int T_index) {
    L = P.Lvec(L_index);
    T = P.Tvec(T_index);
    equilibrium_index = -1;
    max_cycles = P.max_cycles;
    max_chunks = std::round(max_cycles*L*L/target_chunks);
    target_cycles = P.target_cycles;
    target_chunks = std::round(target_cycles*L*L/chunk_sz);
    chunk_sz = P.cunk_sz;
    sample_frequency = P.sample_freqency;
    target_chunks P.target_chunks;
    E = arma::ivec(max_chunks*chunk_sz);
    M = arma::ivec(max_chunks*chunk_sz);
    E_chunks = arma::vec(max_chunks+1);
}

void write_to_file(const Results& r);

class IsingLattice {
    private:
        const double beta // Coldness
        const double J; // Coupling constant of spin
        const int L; // Lattice side length in number of spins
        const arma::vec exp_cache; // Cache of all possible values exp(-beta*DeltaE)
        arma::imat S; // Spin state

        int E; // Total energy
        int M; // Total magnetization

        &int periodic(int k, int l); // Periodic access function

        /* (Re)Calculates energy and magnetization of spin configuration.*/
        void init_E();
        void init_M();
        /* Retrieves cached value of expression exp(-beta*DeltaE)*/
        double exp_DE(int DE);
        /* Advance lattice one cycle. */
        void advance();
    public:
        /* Construct an LxL Ising lattice in an ordered spin configuration with
        all spins up in an environment at temperature T, with coupling constant J. */
        IsingLattice(int L, double T);

        /* Randomize the spins */
        void randomize();

        void run(Results& R, RunParameters RP);
};

IsingLattice::IsingLattice(int L_in, double T_in,)
    : beta(1/T_in),
      L(L_in),
      J(1.0),
      E(0.0),
      M(0.0),
      exp_cache {
        std::exp(-beta*J*8),
        std::exp(-beta*J*4),
        std::exp(beta*J*0),
        std::exp(beta*J*4),
        std::exp(beta*J*8)
      },
      S(L_in, L_in)
{
    init_E();
    init_M();
}

&int IsingLattice::periodic(k,l) {
    return S((k+L)%L, (l+L)%L);
}

void IsingLattice::init_E() {
    E = 0;
    for (int i = 0; i < L; ++i) {
        for (int j = i%2; j < L; i+=2) {
            int s_ij = S(i,j);
            E += -s_ij*(
                  periodic(i-1, j)
                + periodic(i+1, j)
                + periodic(i, j-1)
                + periodic(i, j+1)
            );
        }
    }
}
void IsingLattice::init_M() {
    M = 0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j){
            M += S(i,j);
        }
    }
}

double IsingLattice::exp_DE(int DE) {
    return exp_cache[(DE+8)/4]
}

void IsingLattice::advance() {
    // Pick a random spin:
    int i = arma::randi(arma::distr_param(0,L-1));
    int j = arma::randi(arma::distr_param(0,L-1));
    // Calculate change in energy if spin is flipped:
    int s_ij = S(i,j);
    DE = 2*s_ij*(
          periodic(i-1,j)
        + periodic(i+1,j)
        + periodic(i,j-1)
        + periodic(i,j+1)
    );
    // If flipping lowers total energy: flip, else calculate probability of flip
    // happening regardless in keeping with the Boltzmann distribution.
    if ((DE <= 0) || (exp_DE(DE) >= arma::randu())) {
        S(i,j) *= -1;
        E += DE;
        M += -2*s_ij;
    }
}

void IsingLattice::randomize() {
    S = arma::randi(S.n_rows, S.n_cols, arma::distr_param(0,1)).replace(0,-1);
    init_E();
    init_M();
}

void IsingLattice::run(Results& R) {
    // Take an initial sample
    R.E(0) = E;
    R.M(0) = M;
    int chunk_start = 0;
    int chunk_end = R.chunk_sz; // First chunk has length chunk_sz + 1
    for (int i = 0; i < R.max_chunks; ++i) { // Limit number of cycles
        for (int j = 0; j < R.chunk_sz; ++j) { // Take enough samples for a chunk
            for (int k = 0; k < R.sample_frequency; ++k) { // Make enough flips for a sample.
                advance()
            }
            // Take sample
            R.E(1 + j + i*R.chunk_sz) = E; // The 1 represents the initial sample
            R.M(1 + j + i*R.chunk_sz) = M;
        }
        // Calculate chunk average
        R.E_chunks(i) = arma::mean(R.E.subvec(chunk_start, chunk_end));
        chunk_start = chunk_end + 1;
        chunk_end = chunk_start + R.chunk_sz - 1;
        // Check if we have equilibrium
        if (i > R.target_chunks-2) { // Only check if we have sufficient chunks to end
            int first_chunk = i - R.target_chunks;
            int last_chunk = i;
            sigma = arma::stddev(R.E_chunks.subvec(i-R.target_chunks, i));
            if (sigma < R.epsilon*2*J) { // We have sufficient chunks at equilibrium
                // Set start of equilibrium conditions
                R.equilibrium_index = first_chunk * R.chunk_sz;
                break;
            }
        }
    }
}

void multirun(Parameters P) {
    // Ensure results folder exists
    std::filesystem::create_directory(data_directory);

    n_T = P.Tvec.n_elem;
    n_L = P.Lvec.n_elem;
    n_sets = n_T * n_L;
    n_per_set = P.n_runs_per_set
    n_total_runs =  n_parameter_sets * P.n_per_set;
    std::vector<results> R(n_runs);
    std::cout << "Running...\n";
    for (int i = 0; i < n_L; ++i)
        for (int j = 0; j < n_T; ++j)
            for (int k = 0; k < n_per_set; ++k) {
                lat = IsingLattice(P.Lvec(i), P.Tvec(j));
                if (P.randomize)
                    lat.randomize();
                int l = i*n_T*n_per_set + j*n_per_set + k;
                R[l].initialize(P, i, j);
                lat.run(R[l]);
            }
    std::cout << "Writing...\n";
    for (auto it = R.begin(); it != R.end(); ++it) {
        write_to_file(*it);
    }

}

void write_to_file(const Results& r) {
    // Find an availabel file name
    std::filesystem::path dir(data_directory);
    std::filesystem::path p;
    string base_name = "output_"
    string file_ending = ".txt"
    for (int i = 1; true; ++i) {
        std::stringstream name;
        name << base_name;
        name << std::setw(4) << std::setfill('0') << i;
        name << std::setw(0) << std::setfill(' ') << file_ending;
        p = dir / name.str();
        if (!std::filesystem::exists(p))
            break;
    }

    // Open and prepare file
    std::fstream file(p, std::ios::out | std::ios::trunc);
    int d = std::numeric_limits<double>::digits10; // Use all the precision we have
    output << std::scientific << std::setprecision(d);

    // Write
}
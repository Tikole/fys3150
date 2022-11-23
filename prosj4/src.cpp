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

Parameters::Parameters() {
    runs = 1;
    randomize = true;
    max_chunks = 100;
    target_chunks = 10;
    chunk_sz = 1000;
    sample_frequency = 5;
    epsilon = 0;
}

/* Data structure serving to pass arguments to a run and return results thereof*/
struct Results {
    int L;
    double T;
    bool random;
    int equilibrium_index;
    double max_cycles;
    int max_chunks;
    double target_cycles;
    int target_chunks;
    int chunk_sz;
    int sample_frequency;
    double epsilon;
    arma::ivec E;
    arma::ivec M;
    arma::vec E_chunks;
    arma::ivec chunk_index;
    
    void initialize(const Parameters& P, int T_index);
};

void Results::initialize(const Parameters& P, int T_index) {
    L = P.L;
    T = P.Tvec(T_index);
    random = P.randomize;
    equilibrium_index = -1; // -1 indicates: equilibrium not reached
    target_chunks = P.target_chunks;
    max_chunks = P.max_chunks;
    target_cycles = 0;
    max_cycles = 0;
    chunk_sz = P.chunk_sz;
    sample_frequency = P.sample_frequency;
    epsilon = P.epsilon;
    E = arma::ivec(max_chunks*chunk_sz + 1);
    M = arma::ivec(max_chunks*chunk_sz + 1);
    E_chunks = arma::vec(max_chunks);
    chunk_index = arma::ivec(max_chunks);
}

void write_to_file(const Results& r);

class IsingLattice {
    private:
        const double beta; // Coldness
        const double J; // Coupling constant of spin
        const int L; // Lattice side length in number of spins
        const arma::vec exp_cache; // Cache of all possible values exp(-beta*DeltaE)
        arma::imat S; // Spin state

        int E; // Total energy
        int M; // Total magnetization

        int periodic(int k, int l); // Periodic access function

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

        void run(Results& R);
};

IsingLattice::IsingLattice(int L_in, double T_in)
    : beta(1/T_in),
      L(L_in),
      J(1.0),
      E(0.0),
      M(0.0),
      exp_cache {
        std::exp(beta*J*8),
        std::exp(beta*J*4),
        std::exp(beta*J*0),
        std::exp(-beta*J*4),
        std::exp(-beta*J*8)
      },
      S(L_in, L_in, arma::fill::ones)
{
    init_E();
    init_M();
}

int IsingLattice::periodic(int k, int l) {
    return S((k+L)%L, (l+L)%L);
}

void IsingLattice::init_E() {
    E = 0;
    for (int i = 0; i < L; ++i) {
        for (int j = i%2; j < L; j+=2) {
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
    return exp_cache((DE+8)/4);
}

void IsingLattice::advance() {
    // Pick a random spin:
    int i = arma::randi(arma::distr_param(0,L-1));
    int j = arma::randi(arma::distr_param(0,L-1));
    // Calculate change in energy if spin is flipped:
    int s_ij = S(i,j);
    int DE = 2*s_ij*(
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
                advance();
            }
            // Take sample
            R.E(1 + j + i*R.chunk_sz) = E; // The 1 represents the initial sample
            R.M(1 + j + i*R.chunk_sz) = M;
        }
        // Calculate chunk average
        R.E_chunks(i) = arma::mean(R.E.subvec(chunk_start, chunk_end));
        R.E_chunks(i) /= (L*L); // For purpose of calculation normalize to per spin.
        R.chunk_index(i) = std::round(0.5 * (chunk_start + chunk_end));
        chunk_start = chunk_end + 1;
        chunk_end = chunk_start + R.chunk_sz - 1;
        // Check if we have equilibrium, but only if we have reached target
        // and epsilon has not been set to 0.
        if (i > R.target_chunks-2 && R.epsilon > 0) {
            int first_chunk = i - R.target_chunks + 1;
            int last_chunk = i;
            double variance = arma::var(R.E_chunks.subvec(first_chunk, last_chunk));
            if (variance < R.epsilon*2*J) { // We have sufficient chunks at equilibrium
                // Set start of equilibrium conditions
                R.equilibrium_index = first_chunk * R.chunk_sz;
                // Shed rows
                R.E.shed_rows(2+i*R.chunk_sz, R.E.n_rows-1);
                R.M.shed_rows(2+i*R.chunk_sz, R.M.n_rows-1);
                R.E_chunks.shed_rows(last_chunk+1, R.E_chunks.n_rows-1);
                R.chunk_index.shed_rows(last_chunk+1, R.chunk_index.n_rows-1);
                break;
            }
        }
    }
    R.E_chunks *= (L*L); // Re-multiply with N for consistency with E and M.
}

void multirun(Parameters P) {
    // Check parameters for logic:
    if (P.L < 1 || P.Tvec.n_elem < 1) {
        throw std::invalid_argument("L and T must be provided");
    }

    // Ensure results folder exists
    std::filesystem::create_directory(data_directory);

    int n_T = P.Tvec.n_elem;;
    int runs = P.runs;
    int total_runs =  n_T * runs;
    std::vector<Results> R(total_runs);
    std::cout << "Running...\n";
#ifdef _OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (int i = 0; i < n_T; ++i) {
        for (int j = 0; j < runs; ++j) {
            IsingLattice lat(P.L, P.Tvec(i));
            if (P.randomize)
                lat.randomize();
            int l = i*runs + j;
            R[l].initialize(P, i);
            lat.run(R[l]);
        }
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
    std::string base_name = "output_";
    std::string file_ending = ".txt";
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
    file << std::scientific << std::setprecision(d);

    // Write
    file << "L T random max_cycles target_cycles chunk_sz sample_freq epsilon "
         << "max_chunk target_chunks" << std::endl;
    file << r.L << ' ' << r.T << ' ' << r.random << ' ' << r.max_cycles << ' '
         << r.target_cycles << ' ' << r.chunk_sz << ' ' << r.sample_frequency
         << ' ' << r.epsilon << ' ' << r.max_chunks << ' ' << r.target_chunks
         << std::endl;
    file << "E_chunk_avg \\n E_chunk_index \\n E \\n M" << std::endl;
    for (int i = 0; i < r.E_chunks.n_elem; ++i) {
        file << r.E_chunks(i) << ' ';
    }
    file << std::endl;
    for (int i = 0; i < r.chunk_index.n_elem; ++i) {
        file << r.chunk_index(i) << ' ';
    }
    file << std::endl;
    for (int i = 0; i < r.E.n_elem; ++i) {
        file << r.E(i) << ' ';
    }
    file << std::endl;
    for (int i = 0; i < r.M.n_elem; ++i) {
        file << r.M(i) << ' ';
    }
    file << std::endl;
    file << "equilibrium_index" << std::endl;
    file << r.equilibrium_index << std::endl;
}

// Needed to see content of arma objects in debugger
void pa(const arma::mat& M) {
    M.print();
}
void pa(const arma::vec& v) {
    v.print();
}
void pa(const arma::ivec& v) {
    v.print();
}
void pa(const arma::imat& M) {
    M.print();
}
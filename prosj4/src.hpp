#include <armadillo>

// Structure holding parameters for a call to multirun() 
struct Parameters {
    int L;
    arma::vec Tvec;
    bool randomize;
    int runs;
    double max_chunks;
    double target_chunks;
    int chunk_sz;
    int sample_frequency;
    double epsilon;

    // Constructs default parameters
    Parameters();
};

/*
Runs multiple Markov Chain Monte Carlo calculations on Ising lattices, and
writes to disk. */
void multirun(Parameters P);
/* Parameters members:
     .L : Lattice side length [number of spins]
     .Tvec : Temperatures [J/k_b]
     .ramdomize : true -> random initial state, false -> all-spins-up initial state
     .runs : Indenpendent runs to make for each combination of T and L
     .max_chunks : Limits number of chunks if equilibrium is not reached.
     .target_cycles : Target number of chunks to run after reaching equilibrium for each run
     .chunk_sz : Number of samples to average for determining equilibrium
     .sample_frequency : Number of flips to make for every sample
     .epsilon : Width of epsilon sausage, the system will be deemed equilibrated when
                std. dev. of a certain number of chunks' average energy per spin falls
                within. Units of smallest possible change to energy of lattice. If zero
                each lattice is run for max_cycles.
*/

#include <armadillo>

// Structure holding parameters for a call to multirun() 
struct Parameters {
    arma::ivec Lvec;
    arma::vec Tvec;
    bool randomize;
    int max_cycles;
    int runs_per_set;
    int chunk_sz;
    int target_chunks;
    double epsilon;

    // Constructs default parameters
    Parameters(
        randomize = true;
        max_cycles = 1000;
        runs_per_set = 8;
        chunk_sz = 50;
        target_chunks = 1000;
        epsilon = 1e-4;
    );
}

/*
Runs multiple Markov Chain Monte Carlo calculations on Ising lattices, and
writes to disk. */
void multirun(Parameters P);
/* Parameters members:
     .Lvec : Lattice side lengths [number of spins]
     .Tvec : Temperatures [J/k_b]
     .ramdomize : true -> random initial state, false -> all-spins-up initial state
     .max_cycles : Limits number of cycles if equilibrium is not reached.
     .runs_per_set : Indenpendent runs to make for each combination of T and L
     .chunk_sz : Number of cycles to average for determining convergence
     .target_chunks : Target number of chunks to sample after reaching equilibrium
     .epsilon : width of epsilon sausage, the system will be deemed equilibrated when
                std. dev. of a certain number of chunks' average energy per spin falls
                within. Units of smallest possible change to energy of lattice. If zero
                each lattice is run for exactly max_cycles.

*/

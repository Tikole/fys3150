#pragma once

#include<armadillo>

#include<string>
#include<vector>

const double Ca_mass = 40.078; // [u]
const double Ca_charge = 1.0; // [e]

const double def_B = 1.0; // [T]
const double def_V = 25; // [mV]
const double def_d = 500; // [um]
const double def_x0 = 20.0; // [um]
const double def_y0 = 0.0;
const double def_z0 = 20.0;
const double def_vx0 = 0.0; // [um/us]
const double def_vy0 = 25.0;
const double def_vz0 = 0.0; 

class Particle {
    public:
        double m; // Mass [u]
        double q; // Charge [e]
        arma::vec r; // Position [um]
        arma::vec v; // Velocity [um/us]
        arma::vec a; // Acceleration [um/us^2]

        Particle(double mass, double charge, arma::vec position, arma::vec velocity);
        std::string s() const;
};

typedef std::vector<Particle> Particles;
typedef std::vector<Particle>::iterator particle_it;
typedef std::vector<Particle>::const_iterator const_particle_it;

class PenningTrap;
typedef std::vector<PenningTrap> PT_states;
typedef std::vector<PenningTrap>::iterator state_it;
typedef std::vector<PenningTrap>::const_iterator const_state_it;

/* Represents state of a Penning trap. */
class PenningTrap {
    public:
        bool coulomb; // Include Coulomb force?
        double t; // Time of state [us]
        double B; // Strength of magnetic field [u/(us e)]
        double V; // Electric potential [(u um^2)/(us^2 e)]
        double d; // Characteristic length [um]
        int np; // Number of particles in simulations. 
        Particles p;

        /* Updates next with state of trap after a single step of length h[us]*/
        // Uses this states derivatives. (Euler step)
        void step(double h, PenningTrap& next) const; 
        // Uses state <with>'s derivatives. (Runge-Kutta step)
        void step(double h, PenningTrap& next, const PenningTrap& with) const;
        // Uses explicitly stated derivatives.
        void step(double h, PenningTrap& next, const arma::mat& V, const arma::mat& A) const;

        /* Extracts acceleration/velocity of this states particles and adds/assigns them
        into matrix A/V, optionally multiplied with a constant c. Needed to calculate mean
        derivatives for Runge Kutta method. A better solution would be to store particle velocity
        and acceleration as matrixes in PenningTrap class to begin with, but this would make the
        Particle class obsolete. For addition compatibility of dimensions is assumed. */
        void add_acceleration(arma::mat& A, double c=1.0) const;
        void add_velocity(arma::mat& V, double c=1.0) const;
        void assign_acceleration(arma::mat& A, double c=1.0) const;
        void assign_velocity(arma::mat& V, double c=1.0) const;

        void accel(); // Calculate acceleration of particles.
    public:
        /* Constructors:*/
        // Constructs analytical test case: 
        PenningTrap();
        /* Construct Penning trap with particles.
        - B0: magnetic field strength in T
        - V0: voltage in mV
        - d0: characteristic length in um */
        PenningTrap(double B0, double V0, double d0, Particles, bool coulomb=true);

        /* Iterators for examining particles */
        particle_it begin();
        particle_it end();
        const_particle_it cbegin() const;
        const_particle_it cend() const;

        /* Integrators.
        Return (steps + 1) states of system for times t, t+h, t+2h,..., t+Dt. */
        PT_states advance_forward_euler(double Dt, int steps) const;
        PT_states advance_runge_kutta_4(double Dt, int steps) const;

        /* Examine system: */
        // Examine electromagnetic field.
        arma::vec MF(const arma::vec& r) const; // Magnetic field
        arma::vec MF(double x, double y, double z) const;
        arma::vec EF(const arma::vec& r) const; // Electric field
        arma::vec EF(double x, double y, double z) const;

        bool is_coulomb() const; // Are coulomb interactions in effect?
        int n_particles() const; // Number of particles in simulation
        double time() const; // Time of state

        /* Writes the state to an output filestream. */
        void write_state(std::fstream& stream) const;
};

/* Write <data> to file. <method> should state which integration method was used.
"RK4" or "FE". 
File is given a name PT_<n_particles>_<time[us]>_<steps>_<coulomb>_<method>_<ident>.txt" e.g.
    '100_50_10000_1_RK4.txt'.
'time' is rounded to nearest int.
<ident> should be a short '_' without whitespace or '_'. Prevents unintentional overwriting of other
files with same number of particles, time, steps etc.
<desc> should be a SINGLE human readable line without any newlines that describes trap,
purpose etc.
Format of file is one line description, one line intended for humans, then any number of lines
ofwhitespace separated floating point numbers with each line corresponding to a single
timepoint. Numbers in order:
    t[us] B[T] V[mV] d[um] x1[um] y1[...] z1 vx1[um/us] vy1[...] vz1 x2 y2 z2 vx2 vy2 vz2 ...
 */
void write_to_file(const PT_states& data, std::string method, std::string ident, std::string desc);
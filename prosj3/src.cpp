#include "src.hpp"

#include<cmath>
#include<fstream>
#include<iomanip>
#include<limits>

/* Conversion factors to internal units. */
const double T_conversion_factor = 9.64852558e1; // From tesla
const double mV_conversion_factor = 9.64852558e4; // From millivolt
/*Coulomb constant in [u um^3/us^2 e^2]*/
const double k_e = 1.38935333e5;
/* Pi */
const double pi = 3.14159265358979323846;

// *** class Particle ***
Particle::Particle()
    : m(1), q(1)
{

}
Particle::Particle(double mass, double charge, arma::vec position, arma::vec velocity)
    : m(mass), q(charge), r(position), v(velocity) /* a is initiallized by PenningTrap*/
{}
Particle::Particle(double mass, double charge, double sr, double sv)
    : m(mass), q(charge)
{
    r = arma::vec(3).randn() * sr;
    v = arma::vec(3).randn() * sv;
}
std::string Particle::s() const {
    std::stringstream ss;
    ss << "r: (" << r[0] << ' ' << r[1] << ' ' << r[2] << "), "
        << "v: (" << v[0] << ' ' << v[1] << ' ' << v[2] << "), "
        << "a: (" << a[0] << ' ' << a[1] << ' ' << a[2] << ")";
    return ss.str();
}

// *** class PenningTrap ***
PenningTrap::PenningTrap()
    : t(0.0), coulomb(false), B(def_B*T_conversion_factor), V(def_V*mV_conversion_factor), d(def_d),
      p({Particle(Ca_mass, Ca_charge, {def_x0, def_y0, def_z0}, {def_vx0, def_vy0, def_vz0})}),
      np(1)
{
    p.shrink_to_fit();
    update_accel();
    update_voltage();
}
PenningTrap::PenningTrap(double B0, double V0, double d0, Particles pv, bool b)
    : t(0.0), coulomb(b), B(B0*T_conversion_factor), d(d0), p(pv), np(pv.size()),
      VF_V_mean(V0*mV_conversion_factor), VF_A(0), VF_om(0), VF_ph(0)
{
    p.shrink_to_fit(); // Number of particles is constant.
    update_accel(); // Set initial acceleration of particles.
    update_voltage();
}
PenningTrap::PenningTrap(double B0, double V0, double VF_freq, double VF_Amp,
                         double VF_phase, double d0, Particles pv, bool b)
    : t(0.0), coulomb(b), B(B0*T_conversion_factor), VF_V_mean(V0*mV_conversion_factor), d(d0),
      p(pv), np(pv.size()), VF_A(VF_Amp), VF_om(2*pi*VF_freq),
      VF_ph(VF_phase) 
{
    p.shrink_to_fit(); // Number of particles is constant.
    update_accel(); // Set initial acceleration of particles.
    update_voltage();
}

// Iterators over particles
const_particle_it PenningTrap::cbegin() const {
    return p.cbegin();
}
const_particle_it PenningTrap::cend() const {
    return p.cend();
}

// Integrator functions
PT_states PenningTrap::advance_forward_euler(double Dt, int steps) const {
    int N = steps + 1;
    double h = Dt/steps;
    PT_states V(N, *this); // Sequence of states, our return object
    for (auto it = V.begin() + 1; it != V.end(); ++it) {
        (it-1)->step(h, *it); // Create this state by stepping with previous state.
    }
    return V;
}
PT_states PenningTrap::advance_runge_kutta_4(double Dt, int steps) const {
    int N = steps + 1;
    double h = Dt/steps;
    PT_states V(N, *this); // Sequence of states, our return object
    /* Matrixes for calculating weigted average acceleration and velocity of states.*/
    arma::mat Acc;
    arma::mat Vel;
    for (auto it = V.begin() + 1; it != V.end(); ++it) {
        const PenningTrap& S1 = *(it-1);
        PenningTrap S2; S1.step(h/2, S2);
        PenningTrap S3; S1.step(h/2, S3, S2);
        PenningTrap S4; S1.step(h, S4, S3);
        S1.assign_acceleration(Acc, 1.0/6.0); S1.assign_velocity(Vel, 1.0/6.0);
        S2.add_acceleration(Acc, 2.0/6.0); S2.add_velocity(Vel, 2.0/6.0);
        S3.add_acceleration(Acc, 2.0/6.0); S3.add_velocity(Vel, 2.0/6.0);
        S4.add_acceleration(Acc, 1.0/6.0); S4.add_velocity(Vel, 1.0/6.0);

        (it-1)->step(h, *it, Vel, Acc);
    }
    return V;
}

// Subsidiary functions of integrators
void PenningTrap::step(double h, PenningTrap& next) const {
    step(h, next, *this);
}
void PenningTrap::step(double h, PenningTrap& next, const PenningTrap& with) const {
    arma::mat A; with.assign_acceleration(A);
    arma::mat V; with.assign_velocity(V);
    step(h, next, V, A);
}
void PenningTrap::step(double h, PenningTrap& next, const arma::mat& V, const arma::mat& A) const {
    next = *this;
    next.t = t + h;
    for (int i = 0; i < np; ++i) {
        next.p[i].r = p[i].r + V.col(i)*h;
        next.p[i].v = p[i].v + A.col(i)*h;
    }
    next.update_accel();
    next.update_voltage();
}

void PenningTrap::add_acceleration(arma::mat& A, double c) const {
    for (int i = 0; i < np; ++i) {
        A.col(i) += c*p[i].a;
    }
}
void PenningTrap::add_velocity(arma::mat& V, double c) const {
    for (int i = 0; i < np; ++i) {
        V.col(i) += c*p[i].v;
    }
}
void PenningTrap::assign_acceleration(arma::mat& A, double c) const {
    A.set_size(3, np);
    for (int i =0; i < np; ++i) {
        A.col(i) = c*p[i].a;
    }
}
void PenningTrap::assign_velocity(arma::mat& V, double c) const {
    V.set_size(3, np);
    for (int i = 0; i < np; ++i) {
        V.col(i) = c*p[i].v;
    }
}

// Acceleration:
void PenningTrap::update_accel() {
    for (int i = 0; i < np; ++i) {
        // Initialize force vector with Lorentz force.
        arma::vec F = p[i].q * (EF(p[i].r) + arma::cross(p[i].v, MF(p[i].r)));
        if (coulomb) {
            for (int j = 0; j < np; ++j) {
                if (i != j) {
                    arma::vec R = (p[i].r - p[j].r); // Vector from 1 to 2.
                    double R_sq = arma::dot(R,R);
                    R = arma::normalise(R);
                    F += R*k_e*(p[i].q)*(p[i].q)/R_sq;
                }
            }
        }
        p[i].a = F/(p[i].m);
    }
}
void PenningTrap::update_voltage() {
    V = VF_V_mean * (1 + VF_A*std::cos(VF_om*t + VF_ph));
}

// Examination functions
arma::vec PenningTrap::MF(const arma::vec& r) const{
    return (arma::dot(r,r) < d*d)
        ? (arma::vec{0.0, 0.0, B}) : (arma::vec{0.0, 0.0, 0.0});
}
arma::vec PenningTrap::MF(double x, double y, double z) const {
    return MF({x,y,z});
}
arma::vec PenningTrap::EF(const arma::vec& r) const{
    double f = V/(d*d);
    return (arma::dot(r,r) < d*d)
        ? (arma::vec{f*r(0), f*r(1), -2*f*r(2)}) : (arma::vec{0.0, 0.0, 0.0});
}
arma::vec PenningTrap::EF(double x, double y, double z) const {
    return EF({x,y,z});
}

bool PenningTrap::is_coulomb() const {
    return coulomb;
}
int PenningTrap::n_particles() const {
    return np;
}
double PenningTrap::time() const {
    return t;
}
int PenningTrap::n_trapped() const {
    int n = 0;
    for (auto it = p.cbegin(); it != p.cend(); ++it) {
        n += (arma::dot(it->r, it->r) < d*d);
    }
    return n;
}

// Output functions:
void PenningTrap::write_state(std::fstream& stream) const {
    /* We leave formatting to caller. */
    stream << t << ' ' << B << ' ' << V << ' ' << d;
    for (auto it = p.begin(); it != p.end(); ++it) {
        stream << ' ' << it->r(0) << ' ' << it->r(1) << ' ' << it->r(2) << ' '
               << it->v(0) << ' ' << it->v(1) << ' ' << it->v(2);
    }
    stream << std::endl;
}

// *** Stand alone functions: ***
void write_to_file(const PT_states& data, std::string method, std::string ident, std::string desc) {
    int n_particles = data.front().n_particles();
    int time = std::round(data.back().time());
    int N = data.size() - 1;
    int b = data.front().is_coulomb() ? 1 : 0;

    std::string filename = 
        "PT_" + std::to_string(n_particles)
        + '_' + std::to_string(time)
        + '_' + std::to_string(N)
        + '_' + std::to_string(b)
        + '_' + method + '_' + ident  + ".txt";
    std::fstream file(filename, std::ios::trunc | std::ios::out);
    // Scientific format with all available precision.
    int d = std::numeric_limits<double>::digits10;
    file << std::scientific << std::setprecision(d);
    // Write description:
    file << desc << std::endl;
    // Write line descibing format.
    file << "t[us] B[T] V[mV] x1[um] y1 z1 vx1[um/us] vy1 vz1 x2 y2 z2 vx2 vy2 vz2 ..."
         << std::endl;
    // Write each particle in turn.
    const_state_it it = data.begin();
    while (it != data.end()) {
        (it++)->write_state(file);
    }              
}

// Needed to see content of arma objects in debugger
void pa(const arma::mat& M) {
    M.print();
}
void pa(const arma::vec& v) {
    v.print();
}
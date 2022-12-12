#include "slitex.hpp"

#include <iostream>
#include <cmath>
#include <string>

void print_usage();

int main(int argc, char* argv[]){
    SEParams params;
    double Dt = 2.5e-5; // default timestep
    double T = 0.008; // Default sim time
    int i = 1;
    while (i < argc) {
        std::string s(argv[i++]);
        // help
        if (s == "-help") {
            print_usage();
            return 0;
        }
        // Size of grid, number of points along side
        else if (s == "-M") {
            int n = std::stoi(argv[i++]);
            params.M = n; 
        }
        // Size of grid, distance between points
        else if (s == "-h") {
            double d = std::stod(argv[i++]);
            params.M = std::round(1.0/d) + 1;
        }
        // Timestep as single float
        else if (s == "-Dt") {
            double d = std::stod(argv[i++]);
            Dt = d;
        }
        // Time to simulate as float
        else if (s == "-T") {
            double d = std::stod(argv[i++]);
            T = d;
        }
        // Number of slits to put in wall, 0 also means no wall.
        else if (s == "-slits") {
            int n = std::stoi(argv[i++]);
            params.n_slits = n;
        }
        // Size of slits
        else if (s == "-ap") {
            double d = std::stod(argv[i++]);
            params.aperture = d;
        }
        // Distance between slits
        else if (s == "-sep") {
            double d = std::stod(argv[i++]);
            params.separation_width = d;
        }
        // Initial position of wavepacket 2 floats x and y
        else if (s == "-pos") {
            double d1 = std::stod(argv[i++]);
            double d2 = std::stod(argv[i++]);
            params.x0 = d1;
            params.y0 = d2;
        }
        // Initial width of wavepacket (std.dev), 2 floats x and y
        else if (s == "-sigma") {
            double d1 = std::stod(argv[i++]);
            double d2 = std::stod(argv[i++]);
            params.sigma_x = d1;
            params.sigma_y = d2;
        }
        // Initial momentum of wavepacket, 2 floats px, py
        else if (s == "-mom") {
            double d1 = std::stod(argv[i++]);
            double d2 = std::stod(argv[i++]);
            params.px = d1;
            params.py = d2;
        }
    }
    run(params, Dt, T);

    return 0;
}

void print_usage() {
    std::cout << "Run an electron diffraction experiment.\nResults outputted to "
                 "a file 'slitex.hdf5'. Following arguments accepted:\n"
                 "    -help       - prints this message.\n"
                 "    -M          - Number of points along grid. Default: 200.\n"
                 "    -h          - Distance between points of grid. Default 1/200.\n"
                 "                  Alternative way of setting M. Do not set both.\n"
                 "    -Dt         - Length of time step. Default: 0.000025.\n"
                 "    -T          - Length of time to simulate. Default: 0.008.\n"
                 "    -slits      - Number of slits to make, 0 also means no wall.\n"
                 "                  Default: 0.\n"
                 "    -ap         - Size of slits. Default: 0.05.\n"
                 "    -sep        - Width of wall between slits. Default: 0.05.\n"
                 "    -pos        - Initial position of center of wave packet.\n"
                 "                  Two arguments x and y. Default: 0.25 0.5\n"
                 "    -sigma      - Spread in initial wavepacket, std.dev. Two \n"
                 "                  arguments. Default 0.05 0.05\n"
                 "    -mom        - Initial momentum of particle. 2 arguments.\n"
                 "                  Default: 200.0 0.0\n";
}
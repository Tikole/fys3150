#pragma once

struct SEParams {
    int M; // Sidelength of grid in number of points. Barrier included.
    int n_slits; // Number of slits 0 (also disables wall), 1, 2, 3.
    double wall_thickness; // Thickness of center wall.
    double wall_position; // Position of center wall along x-axis.
    double separation_width; // Width of wall pieces separating slits.
    double aperture; // Aperture of slits.
    double V0; // Potential inside barriers. 
    double x0; // Initial x-position of center of wave packet.
    double y0; // Initial y-position ...
    double sigma_x; // std. dev. of gaussian envelope of initial wave packet along x-axis.
    double sigma_y; // ... along y-axis.
    double px; // Initial momentum of wave packet along x-axis.
    double py; // ... along y-axis.

    /* Constructs defaults*/
    SEParams();
};

/* Runs a electron interference experiment and writes results to disk as 'slitex.hdf5'*/
void run(const SEParams& params, double Dt, double T);
#include "slitex_impl.hpp"

/* Make the A and B matrix defining the time evolution of the wave packet */
void make_A_and_B(const mat& V, double h, double Dt, sp_cx_mat& A_out, sp_cx_mat& B_out) {    
    assert(V.is_square());
    assert(h > 0.0);
    assert(Dt > 0.0);
    int n = V.n_cols;
    int N = V.n_elem;
    cx_double r(0, Dt/(2*h*h));
    cx_double s(0, Dt/2);
    A_out.zeros(N, N);
    A_out.diag() = 1.0 + 4.0*r + s*V.as_col();
    A_out.diag(-1).fill(-r);
    A_out.diag(1).fill(-r);
    for (int i=n-1; i<A_out.diag(1).n_elem; i+=n) {
        A_out.diag(-1)(i) = cx_double(0.0, 0.0);
        A_out.diag(1)(i) = cx_double(0.0, 0.0);
    }
    A_out.diag(-n).fill(-r);
    A_out.diag(n).fill(-r);

    B_out.zeros(N,N);
    B_out.diag() = 1.0 - 4.0*r - s*V.as_col();
    B_out.diag(-1).fill(r);
    B_out.diag(1).fill(r);
    for (int i=n-1; i<B_out.diag(1).n_elem; i+=n) {
        B_out.diag(-1)(i) = cx_double(0.0, 0.0);
        B_out.diag(1)(i) = cx_double(0.0, 0.0);
    }
    B_out.diag(-n).fill(r);
    B_out.diag(n).fill(r);
}
/* Set up the matrix discretization of the potential*/
void make_V(const SEParams& params, mat& V_out) {
    int M = params.M;
    int n_slits = params.n_slits;
    const double& V0 = params.V0;
    V_out.zeros(M, M);
    // If we are to have slits set up central partition with slits
    if (n_slits > 0) {
        const double& pos = params.wall_position;
        const double& thick = params.wall_thickness;
        const double& ap = params.aperture;
        const double& sep = params.separation_width;
        double h = 1.0/(M-1);
        // Finding wall limits
        double x0 = pos - 1.0/2*thick; // Start of wall
        double x1 = pos + 1.0/2*thick; // End of wall
        // Column index range of wall
        int row0 = std::round(x0/h);
        int row1 = std::round(x1/h);
        // Set the whole wall to V0, then carve out slits.
        V_out.rows(row0, row1).fill(V0);
        double y0; // y of first cut
        if (n_slits%2) { // Odd number of slits
            y0 = 0.5*(1-ap);
            ++n_slits; // Since the first cut in the odd case is done twice this
                       // lets us treat odd and even cuts the same.
        }
        else { // Even number of slits
            y0 = 0.5*(1+sep);
        }
        for (int i=0; 2*i < n_slits; ++i) {
            double ya = y0 + i*(ap+sep);
            double yb = ya + ap;
            int col_a = std::round(ya/h);
            int col_b = std::round(yb/h);
            V_out.cols(col_a, col_b).fill(0.0);
            // Flip indexes around center and make another cut
            col_a = V_out.n_cols - 1 - col_a;
            col_b = V_out.n_cols - 1 - col_b;
            V_out.cols(col_b, col_a).fill(0.0);
        }
    }
    // Set up outer walls
    V_out.row(0).fill(V0);
    V_out.row(M-1).fill(V0);
    V_out.col(0).fill(V0);
    V_out.col(M-1).fill(V0);
}
/* Set up the initial state of the wave function, called S0 here rather than U0. */
void make_S0(const SEParams& params, cx_vec& S_out) {
    const int& M = params.M;
    const double h = 1.0/(M-1);
    const double& xc = params.x0;
    const double& yc = params.y0;
    const double& sx = params.sigma_x;
    const double& sy = params.sigma_y;
    const double& px = params.px;
    const double& py = params.py;
    cx_mat S;
    S.zeros(M,M);
    // Only set internal elements, edge elements should be exactly zero.
    for (int i = 1; i < S.n_rows-1; ++i)
        for (int j=1; j < S.n_cols-1; ++j) {   
                double x = i*h;
                double y = j*h;
                double re = - ((x-xc)*(x-xc))/(2*sx*sx) - ((y-yc)*(y-yc))/(2*sy*sy);
                double im = px*(x-xc) + py*(y-yc);
                S(i,j) = std::exp(cx_double(re,im));
            }
    S_out = arma::normalise(S.as_col());
}

SEParams::SEParams(){
    M = 200;
    n_slits = 0;
    wall_thickness = 0.02;
    wall_position = 0.5;
    separation_width = 0.05;
    aperture = 0.05;
    V0 = 1.0e10;
    x0 = 0.25;
    y0 = 0.5;
    px = 200.0;
    py = 0.0;
    sigma_x = 0.05;
    sigma_y = 0.05;
}

void run(const SEParams& params, double Dt, double T) {
    auto t0 = std::chrono::high_resolution_clock::now();
    std::cout << "Setting up... ";
    std::cout.flush();
    // Set up V matrix
    mat V;
    make_V(params, V);
    // Set up A and B
    double h = 1.0/(params.M-1);
    sp_cx_mat A;
    sp_cx_mat B;
    make_A_and_B(V, h, Dt, A, B);
    // Set up initial state
    cx_vec S0;
    make_S0(params, S0);
    // States vector
    int n_steps = std::round(T/Dt);
    cx_mat S(S0.n_elem, n_steps+1);
    S.col(0) = S0;
    auto t1 = std::chrono::high_resolution_clock::now();
    double setup_time = std::chrono::duration<double>(t1 - t0).count();
    std::cout << '(' << setup_time << " s)\n";
    std::cout << "Calculating...\n";
    // Solve for S1, S2, ... and write to disk.
    std::cout << '0' << '/' << S.n_cols;
    std::cout.flush();
    for (int i=1; i<S.n_cols; ++i) {
        S.col(i) = arma::spsolve(A, B*S.col(i-1));
        std::cout << "\r          \r";
        std::cout << i << '/' << S.n_cols;
        std::cout.flush();
    }
    std::cout << "\r          \r" << S.n_cols << '/' << S.n_cols << '\n';
    std::cout.flush();
    auto t2 = std::chrono::high_resolution_clock::now();
    double solving_time = std::chrono::duration<double>(t2 - t1).count();
    // Write to disk
    std::cout << " (" << solving_time << "s)\n";
    std::cout << "Writing... ";
    std::cout.flush();
    arma::mat S_real = arma::real(S);
    arma::mat S_imag = arma::imag(S);
    S_real
        .save(arma::hdf5_name("slitex.hdf5", "S_real"));
    S_imag
        .save(arma::hdf5_name("slitex.hdf5", "S_imag", arma::hdf5_opts::append));
    V
        .save(arma::hdf5_name("slitex.hdf5", "V", arma::hdf5_opts::append));
    vec{Dt}
        .save(arma::hdf5_name("slitex.hdf5", "Dt", arma::hdf5_opts::append));
    vec{T}
        .save(arma::hdf5_name("slitex.hdf5", "T", arma::hdf5_opts::append));
    ivec{n_steps+1}
        .save(arma::hdf5_name("slitex.hdf5", "N", arma::hdf5_opts::append));
    ivec{params.M}
        .save(arma::hdf5_name("slitex.hdf5", "M", arma::hdf5_opts::append));
    ivec{params.n_slits}
        .save(arma::hdf5_name("slitex.hdf5", "n_slits", arma::hdf5_opts::append));
    vec{params.wall_thickness}
        .save(arma::hdf5_name("slitex.hdf5", "wall_thickness", arma::hdf5_opts::append));
    vec{params.wall_position}
        .save(arma::hdf5_name("slitex.hdf5", "wall_position", arma::hdf5_opts::append));
    vec{params.separation_width}
        .save(arma::hdf5_name("slitex.hdf5", "separation_width", arma::hdf5_opts::append));
    vec{params.aperture}
        .save(arma::hdf5_name("slitex.hdf5", "aperture", arma::hdf5_opts::append));
    vec{params.V0}
        .save(arma::hdf5_name("slitex.hdf5", "V0", arma::hdf5_opts::append));
    vec{params.x0, params.y0}
        .save(arma::hdf5_name("slitex.hdf5", "pos", arma::hdf5_opts::append));
    vec{params.sigma_x, params.sigma_y}
        .save(arma::hdf5_name("slitex.hdf5", "sigma", arma::hdf5_opts::append));
    vec{params.px, params.py}
        .save(arma::hdf5_name("slitex.hdf5", "p", arma::hdf5_opts::append));
    auto t3 = std::chrono::high_resolution_clock::now();
    double write_time = std::chrono::duration<double>(t3 - t2).count();
    double total_time = std::chrono::duration<double>(t3-t0).count();
    std::cout << '(' << write_time << " s)\n";
    std::cout << "Total time: " << total_time << " s\n";
}

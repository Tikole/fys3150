#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<limits>

double u(double x){
    return 1.0 - (1 - std::exp(-10.0))*x - std::exp(-10.0*x);
}

int main(int argc, char* argv[]){
    /* Number of points to evaluate is a provided as a command line argument.
       If not provided use a default number. */
    const int def_evaluation_points = 1000;
    int n = (argc > 1) ? std::stoi(argv[1]) : def_evaluation_points;
    
    /* Evaluate function between, and including, x0 and x1. */
    const double x0 = 0.0;
    const double x1 = 1.0;
    double dx = (x1 - x0)/double(n-1);
    std::vector<double> X(n);
    int i = 0;
    while (i < X.size()-1) {
        X[i] = x0 + i*dx;
        ++i;
    }
    X[i] = x1;

    std::vector<double> Y(n);
    for (auto i = 0; i < X.size(); ++i) {
        Y[i] = u(X[i]);
    }

    /* Write numbers with all the precision allowed by the double type in scientific format. */
    std::fstream output_file("output.txt", std::ios::trunc | std::ios::out);
    output_file << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
    for (int i = 0; i < X.size(); ++i){
        output_file << X[i] << " " << Y[i] << std::endl;
    }

    return 0;
}
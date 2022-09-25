#include "src.hpp"

#include<iostream>
#include<string>
#include<vector>


void print_usage() {
    std::cout << "NOT YET IMPLEMENTED" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc == 1) { // Only one arg: print usage and quit
        print_usage();
        return 0;
    }
    else{
        std::string command = std::string(argv[1]);
        if (command == std::string("test")) { // RUn all tests, print error counts
            std::cout << "Running tests\n---\n"
                      << "Test 1...\n"
                      << test_1() << " errors.\n";
            std::cout << "Test 2...\n"
                      << test_2() << " errors.\n";
            std::cout << "Test 3...\n"
                      << test_3() << " errors.\n";

        }
        else if (command == std::string("scaling")) { 
            /* Test scaling of our implementation of the jacobi method.
            After the 'scaling' command there should be a single integer
            argument n > 2. Jacobi method is used on matrixes of dimensions
            3x3, 4,4,..., NxN.*/
            if (argc != 3) {
                std::cout << "The 'scaling' command should be followed by a single integer argument.";
            }
            else {
                int n;
                try {
                    n = std::stoi(argv[2]);
                }
                catch(...) {
                    std::cout << "The 'scaling' command should be followed by a single integer argument.";
                    return 0;
                }
                scaling_test(n, "scaling");
                scaling_test(n, "densescaling");
            }
        }
        else {
            std::cout << "Unrecognized command: '" << std::string(argv[1]) << "'. For help provide no arguments.\n";
        }
    }
}
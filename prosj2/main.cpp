#include "src.hpp"

#include<iostream>
#include<string>


void print_usage() {
    std::cout << "NOT YET IMPLEMENTED" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc == 1) { // Only one arg: print usage and quit
        print_usage();
        return 0;
    }
    else{
        if (std::string(argv[1]) == "test") { // RUn all tests, print error counts
            std::cout << "Running tests\n---\n"
                      << "Test 1...\n"
                      << test_1() << " errors.\n";
            std::cout << "Test 2...\n"
                      << test_2() << " errors.\n";

        }
    }
}
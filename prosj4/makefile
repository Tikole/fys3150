all:
	g++ main.cpp src.cpp -std=c++17 -larmadillo -g -O0 -o isingdb
	g++ main.cpp src.cpp -std=c++17 -larmadillo -O3 -o isingST
	g++ main.cpp src.cpp -std=c++17 -larmadillo -fopenmp -O3 -o ising

debug:
	g++ main.cpp src.cpp -std=c++17 -larmadillo -g -O0 -o isingdb
scrub:
	rm data/*.txt
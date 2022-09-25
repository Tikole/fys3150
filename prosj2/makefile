all:
	g++ main.cpp src.cpp -larmadillo -g -O0 -o main.exe
	make test
build:
	g++ main.cpp src.cpp -larmadillo -g -O0 -o main.exe
test:
	./main.exe test
scrub:
	rm *.exe* *scaling*.txt *.pdf
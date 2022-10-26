all:
	make debug
	make final
debug:
	g++ main.cpp src.cpp -larmadillo -g -O0 -o dbmain.exe
final:
	g++ main.cpp src.cpp -larmadillo -O3 -o main.exe
scrub:
	rm a.out PT*txt *.exe *.png *.pdf
clean:
	rm PT*.txt *.png *.pdf

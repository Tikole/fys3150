exe: # This assumes hdf5 headers and archives are located as they are for me. Regular Ubuntu running in WSL. 
	g++ main.cpp slitex.cpp -I/usr/include/hdf5/serial/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -larmadillo -lhdf5 -g -O0 -o slitex
junk:
	g++ junk.cpp slitex.cpp -I/usr/include/hdf5/serial/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -larmadillo -lhdf5 -g -O0 -o junk
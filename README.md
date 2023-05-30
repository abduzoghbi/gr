
## GR: a ray tracing code for simulating light rays around a black hole.


### To Compile:
- Make sure you have gsl and hdf5 libraries.
- Define where those libraries and header files are:
INC=""
LIB=""

- run this to compile: 

g++ -fopenmp -I$INC -O0 -g3 -Wall -c -o *cpp
g++ -fopenmp -I$INC -O0 -g3 -Wall -c *cpp

- link:

g++ -fopenmp -L$LIB -o gr *o -lgsl -lgslcblas -lhdf5 -lhdf5_cpp


Run the code for help

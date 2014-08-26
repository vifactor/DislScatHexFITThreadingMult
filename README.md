Program to calculate vertical threading dislocation density in epitaxial layer
with (0001)-type interface from double-crystal scans.
Simultaneous fit to multiple data files is possible.

# Requirements

### Following programs should be present on your computer
- gcc compiler
- cmake

### Following libraries should be present on your computer
- [Gnu Scientific Library (GSL)](http://www.gnu.org/software/gsl/)
- [Basic Linear Algebra Subprograms (BLAS)](http://www.netlib.org/blas/)

### Following libraries will be installed by cmake
- [Port3](http://www.netlib.org/port/)
- [libconfig](http://www.hyperrealm.com/libconfig/)

# Build
From command line in the program folder type following:

1. mkdir build
2. cd build/
3. cmake ..
4. make

You find executable file **DislScatHexFITThreadingMult** in the *build* folder

## Tested on
Linux Mint 17

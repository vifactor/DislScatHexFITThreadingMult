Program to calculate vertical threading dislocation density in hexagonal epitaxial layer
with (0001)-interfaces from double-crystal scans.
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
5. (optional) ln -s /home/user/folder-with-program/DislScatHexFITThreadingMult/build/DislScatCubFITThreadingMult /home/user/bin/DislScatHexFITThreadingMult

After 4-th step you may find *DislScatHexFITThreadingMult* in the *build/* directory.

After 5-th step, you will be able then run your program from any location simply typing 
*DislScatHexFITThreadingMult* in a terminal window.

# Usage
####Suggested folder structure:
.  
├── config  
│   ├── data.cfg  
│   ├── fit.cfg  
│   ├── sample.cfg  
├── data1.dat  
├── data2.dat  
├── ...  
├── datan.dat  

Data files should mandatory contain (at least) two columns entitled in the header as follows:  
    #   [omega]  [intensity]  
or  
    #   [domega]  [intensity]  
if 0-point is in the peak center. [omega]-column values should be in degrees.

Each added datafile should possess the following record in the *data.cfg* configuration file:  
    {  
        file = "relative_path/data.dat";  
        Q = [ H, K, I, L ];  
        resolution :  
        {  
            x = double_value;  
            z = double_value;  
        };
        I0 = double_value;  
        Ibg = double_value;  
    }  
Do not forget to add record  
    {  
		name = "Data.[index].I0";  
		minVal = double_value;  
		maxVal = double_value;  
	}  
to *fit.cfg* upon addition of a new datafile. 
This will add a new intensity-scale fit variable. 
Normally, you always want to fit intensity scale.

####To run the program type:
*DislScatHexFITThreadingMult* config/

####Folder structure after run:  
.
├── config  
│   ├── data.cfg  
│   ├── data.~cfg  
│   ├── fit.cfg  
│   ├── resume.txt  
│   ├── sample.cfg  
│   ├── sample.~cfg  
│   ├── data1.ft  
│   ├── data2.ft  
│   ├── ...  
│   └── datan.ft  
├── data1.dat  
├── data2.dat  
├── ...  
└── datan.dat  

## Tested on
Linux Mint 17
Ubuntu 12.04, 14.04

LCT-INT for spinless fermions
=============================
This the code used for the publication
- Lei Wang, Mauro Iazzi, Philippe Corboz, and Matthias Troyer, "Efficient continuous-time quantum Monte Carlo method for the ground state of correlated fermions", [Phys. Rev. B 91, 235151 (2015)](http://dx.doi.org/10.1103/PhysRevB.91.235151)

It Works for
- Ground state 
- Bipartite lattices 
- Repulsive interaction 

## Prerequisites
- [CMake](https://cmake.org)
- [ALPS](http://alps.comp-phys.org)
- [Eigen3](http://eigen.tuxfamily.org)

## Compilation
    First, create and edit your machine config file in ./src/config/mymachine.cmake, then 
    mkdir build 
    cd build 
    cmake -DUSE_MACHINE=mymachine -DCMAKE_INSTALL_PREFIX=../ -DCMAKE_BUILD_TYPE=Release ../src
    make 
    make install 
This will generate an excutable `../bin/main`.  

## To run
    mkdir ../data 
    mpirun -np 4 ../bin/main  -a 10 -T 120 ../input/params.in 
The result and checkpoint files will be in `../data/`. They are in the [hdf5](https://www.hdfgroup.org/HDF5/) format. 
To inspect them, run 
    h5ls -r ../data/test.out.h5
    h5ls -r ../data/test.chkp/*.h5
To parse and visualize them, you can use the scripts provided in `../analysis`.

**Vola, have fun!**


## Extensions 
- Hubbard model, [Phys. Rev. X 5, 031007 (2015)](http://dx.doi.org/10.1103/PhysRevX.5.031007)
- Mass-imbalanced Hubbard model, [Phys. Rev. B 92, 235129 (2015)](http://dx.doi.org/10.1103/PhysRevB.92.235129)
- Stochastic series expansion, [Phys. Rev. B 93, 155117 (2016)](http://dx.doi.org/10.1103/PhysRevB.93.155117)

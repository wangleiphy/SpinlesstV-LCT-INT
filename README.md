LCT-INT for spinless fermions
=============================
This the code used for the publication 
 - Lei Wang, Mauro Iazzi, Philippe Corboz, and Matthias Troyer, "Efficient continuous-time quantum Monte Carlo method for the ground state of correlated fermions" [Phys. Rev. B 91, 235151 (2015)](http://dx.doi.org/10.1103/PhysRevB.91.235151)

## Prerequisites
- [CMake](https://cmake.org) 
- [ALPS](http://alps.comp-phys.org)
- [Eigen3](http://eigen.tuxfamily.org)

## Compilation
    First, create and edit your machine config file in ./src/config/mymachine.cmake, then 
    mkdir build 
    cd build 
    cmake -DUSE_MACHINE=mymachine -DCMAKE_INSTALL_PREFIX=/install/path/ -DCMAKE_BUILD_TYPE=Release ../src
    make 

## To run
    mkdir ../data 
    mpirun -np 4 ./main  -a 10 -T 120 ../input/params.in 
The results and checkpoint files will be in ../data/ 

## Author 
- Lei Wang (ETH Zurich 2014)

Vola!

LCT-INT (**L**inear scaling **C**ontinous-**T**ime **INT**eraction expansion) for spinless fermions
=============================
This is the (original and fresh) code used for the paper
- Lei Wang, Mauro Iazzi, Philippe Corboz, and Matthias Troyer, "Efficient continuous-time quantum Monte Carlo method for the ground state of correlated fermions", [Phys. Rev. B 91, 235151 (2015)](http://dx.doi.org/10.1103/PhysRevB.91.235151)

It works for
- Ground state 
- Bipartite lattices 
- Repulsive interaction 

## Prerequisites
- [CMake](https://cmake.org)
- [ALPS](http://alps.comp-phys.org)
- [Eigen3](http://eigen.tuxfamily.org)

## Compilation
    First, create and edit your machine specific config file in ./src/config/mymachine.cmake, then 
    mkdir build 
    cd build 
    cmake -DUSE_MACHINE=mymachine -DCMAKE_INSTALL_PREFIX=../ -DCMAKE_BUILD_TYPE=Release ../src
    make 
    make install 
This will generate an excutable `../bin/main`.  

## To run
    mkdir ../data 
    mpirun -np 4 ../bin/main  -a 10 -T 120 ../input/params.in 
The results and checkpoint files will be stored in `../data/`. They are in the [hdf5](https://www.hdfgroup.org/HDF5/) format. 
To inspect them, run 

    h5ls -r ../data/test.out.h5
    h5ls -r ../data/test.chkp/*.h5
To parse and visualize them, you can use the scripts provided in `../analysis/`.

**Vola, have fun!**


## Extensions 
- Hubbard model, [Phys. Rev. X 5, 031007 (2015)](http://dx.doi.org/10.1103/PhysRevX.5.031007)
- Mass-imbalanced Hubbard model, [Phys. Rev. B 92, 235129 (2015)](http://dx.doi.org/10.1103/PhysRevB.92.235129)
- Stochastic series expansion, [Phys. Rev. B 93, 155117 (2016)](http://dx.doi.org/10.1103/PhysRevB.93.155117)

## Run with [Singularity](https://www.sylabs.io/singularity/)
We also provide options to run our program on Singularity. In this way, one should not worry about the difficulty of compiling the code. Only singularity is needed. 

### Step 1:
Build the singularity image (Root or [fakeroot](https://docs.sylabs.io/guides/3.7/user-guide/fakeroot.html) is needed).
```bash
make # You will see main.sif in the current directory
```

### Step 2 (optional):
Modify the input of the program in: `input/params.in` and `input/mylattices.xml` to your need.

### Step 3:
Run the code
```bash
make run # This command will run.
# But if you modify something in Step.2. You may also need to do something in the `Makefile`
```

Especially, the `filename` value in `input/params.in` may be under `$PWD` (or one may need to mount the directory of `filename` in singularity container manually). 

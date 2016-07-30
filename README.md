Code used for the paper http://dx.doi.org/10.1103/PhysRevB.91.235151
Requires:
- CMake https://cmake.org
- ALPS http://alps.comp-phys.org/mediawiki/
- Eigen3 http://eigen.tuxfamily.org

To build:
create and edit your machine config file in ./src/config/mymachine.cmake
mkdir build 
cd build 
cmake -DUSE_MACHINE=mymachine -DCMAKE_INSTALL_PREFIX=/path/to/install/directory -DCMAKE_BUILD_TYPE=Release ../src
make 

To run:
mkdir ../data 
mpirun -np 4 ./main  -a 10 -T 120 ../input/params.in 

Vola!

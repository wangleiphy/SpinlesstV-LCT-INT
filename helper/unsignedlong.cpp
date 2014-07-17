#include <iostream>
int main(){

    unsigned long blocksize=4294967296L; 
    unsigned long itau = 3421126067L; 
    unsigned long b= itau/blocksize; 

    std::cout << b  << std::endl; 
    return 0;
}

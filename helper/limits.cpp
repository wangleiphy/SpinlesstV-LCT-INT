#include <iostream>     // std::cout
#include <limits>       // std::numeric_limits

int main () {

  typedef unsigned T; 

  std::cout << std::boolalpha;
  std::cout << "Minimum value for int: " << std::numeric_limits<T>::min() << '\n';
  std::cout << "Maximum value for int: " << std::numeric_limits<T>::max() << '\n';
  std::cout << "int is signed: " << std::numeric_limits<T>::is_signed << '\n';
  std::cout << "Non-sign bits in int: " << std::numeric_limits<T>::digits << '\n';
  std::cout << "int has infinity: " << std::numeric_limits<T>::has_infinity << '\n';
  return 0;
}

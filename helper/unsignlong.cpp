#include <iostream>     // std::cout
#include <limits>       // std::numeric_limits

int main () {
  std::cout << std::boolalpha;
  std::cout << "Minimum value for int: " << std::numeric_limits<unsigned long>::min() << '\n';
  std::cout << "Maximum value for int: " << std::numeric_limits<unsigned long>::max() << '\n';
  std::cout << "int is signed: " << std::numeric_limits<unsigned long>::is_signed << '\n';
  std::cout << "Non-sign bits in int: " << std::numeric_limits<unsigned long>::digits << '\n';
  std::cout << "int has infinity: " << std::numeric_limits<unsigned long>::has_infinity << '\n';
  return 0;
}

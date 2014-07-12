#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
 
int main()
{
    std::vector<double> data = { 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
 
    auto lower = std::lower_bound(data.begin(), data.end(), 0.25);
    auto upper = std::upper_bound(data.begin(), data.end(), 0.75);
 
    std::copy(lower, upper, std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;  
}

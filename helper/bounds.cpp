#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
 
int main()
{
    std::vector<double> data = {0., 0.2, 1.3 , 2.3, 3.3, 4.3, 5., 6.3, 7.3, 8.1, 9.5, 10};
 
    //auto lower = std::lower_bound(data.begin(), data.end(), 2, std::less_equal<int>()); 
    auto lower = std::lower_bound(data.begin(), data.end(), 0); 
    auto upper = std::upper_bound(data.begin(), data.end(), 5, std::less_equal<double>());
 
    std::copy(lower, upper, std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;  
}

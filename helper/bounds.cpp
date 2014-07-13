#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
 
int main()
{
    std::vector<int> data = {1 , 2, 3, 4, 5, 6, 7, 8, 9};
 
    //auto lower = std::lower_bound(data.begin(), data.end(), 2, std::less_equal<int>()); 3 4 5 6 
    auto lower = std::lower_bound(data.begin(), data.end(), 2); //2 3 4 5 6 

    auto upper = std::upper_bound(data.begin(), data.end(), 6);
 
    std::copy(lower, upper, std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;  
}

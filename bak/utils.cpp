double npickk( unsigned int n, unsigned int k ) // n!/(n-k)! = n*(n-1) * (n-2) * (n-k+1)
{
    if (k > n) return 0.;
    if (k == 0) return 1.;

    double result = n;
    for(unsigned int i = 1; i < k; ++i ) {
        result *= double(n-i);
    }
    return result;
}

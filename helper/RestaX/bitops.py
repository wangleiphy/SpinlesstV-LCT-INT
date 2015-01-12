# written by Matthias Troyer <matthiastroyer@mac.com> for Microsoft Research

def popcount(i):
    '''
    i the 32-bit word of which 1-bits should be counted
    return the number of 1-bits in the word
    '''
    return bin(i).count('1')

def btest(i,p):
    '''
    return true if the p-th bit of i is set, false otherwise
    '''
    return i&(1<<p)

def bset(i,p):
    '''
    return the integer i with the bit at position p set to one
    '''
    return i | (1<<p)

def bget(i,p):
    '''
    return the p-th bit of the word i
    '''
    return (i >> p) & 1

def bclear(i,p):
    '''
    return the integer i with the bit at position p set to zero
    '''
    return i & (~(1<<p))

def bflip(i,p):
    '''
    return the integer i with the bit at position p flipped: (1->0, 0->1)
    '''
    return i ^(1<<p)

def fermionParity(i,p):
  return popcount(i&((1<<p)-1))%2


if __name__=='__main__':
    i = 8
    print bin(i)
    print popcount(i)
    p = 2 
    print bget(i, 2)
    print btest(i, p), (1<<p)
    print bset(i, p), bin(bset(i,p))
    print bclear(i, p)
    print bin(i), p, bin(bflip(i,p))

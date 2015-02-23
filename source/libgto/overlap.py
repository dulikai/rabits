#! /usr/bin/env python
import math

# overlap.py

class overlap:
    """ calc. overlap int. between primary GTOs """
    def __init__(self):
    
        return
        
    @staticmethod    
    def ssu(dist, za, zb):
        """
        calculate un-normalized <s|s> integral
        dist is |A-B|
        zeta are the orbital exponent of pos. A & B;
        """
        zeta = za + zb
        xi = za * zb / zeta
        k = math.pow(math.pi/zeta, 3.0/2.0)
        ek = math.exp(-xi * dist * dist)
        s = k * ek
        return s

        
        
        
if __name__ == "__main__":
    x = overlap()
    for i in xrange(200):
        dist = 1.0 + 0.1*i
        s = x.ssu(dist, 0.1, 0.1)
        print dist, s
    
    
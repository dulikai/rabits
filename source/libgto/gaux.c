#include <stdio.h>
#include <math.h>


/*
    some basic routine for GTO integrals.
    OS scheme.
    see: S. Obara and A. Saika, 
    J. Chem. Phys. 84, 3963 (1986); doi: 10.1063/1.450106
*/

/*
DEFINE MACROs
*/

#ifndef bM_PI
#define bM_PI 3.141592653589793238
#endif  /* bM_PI */

#ifndef bEPS_IN_Fm
#define bEPS_IN_Fm 1.0e-15
#endif  /* bEPS_IN_Fm */

#define PI bM_PI


/*
 overlap integral
 unnormalized Cartesian Gaussian functions
 (a||b) = \int{dr \psi(r;\zeta_a,a,A) \psi(r;zeta_b,b,B)}
 %
 after integral of (s||s) = (\pi/\zeta)^(3/2) exp{-\xi(A-B)^2}
 where \zeta = \zeta_a + \zeta_b & \xi = \zeta_a * \zeta_b / \zeta
 and A B is the nuclear pos. so, we use R = |A-B|.

 <s|s>
*/

// un-normalized <s||s> integral
// $\zeta$ is the orbital exponent
// input: exponent of A,B. and the distance of A,B.
 double 
 overlap_ssu(double dist, double za, double zb)
 {
    double zeta, xi;
    double k, ek, s;
    
    zeta = za + zb;
    xi = za * zb / zeta;
    k = pow(PI/zeta, 3.0/2.0);
    ek = exp(-xi * dist * dist);
    s = k * ek;
    
    return s;
 }
 
 
double 
ovlpu(double R, double a, double b)
{
   double c, ck, zeta, xi;
   double s;

   zeta = a + b;
   xi = a * b/ zeta;
  
   c = pow(PI/zeta, 3.0/2.0);
   ck = exp(-xi * R * R);
   s = c * ck;
   return s;
}


#undef PI

int main()
{
    double R = 1.0;
    double za = 0.1, zb = 0.1;
    double x;
    int i;
    for (i = 0; i < 100; i++) {
        x = overlap_ssu(R+0.1*i, za, zb);
        y = ovlpu(R+0.1*i, za, zb);
        printf("%12.8g%12.8g\n", x, y);
    }
    return 0;
}

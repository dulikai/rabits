#include <stdio.h>
#include <math.h>

#include "gaux.h"

/*
    some basic routine for GTO integrals.
    OS scheme.
    see: S. Obara and A. Saika, 
    J. Chem. Phys. 84, 3963 (1986); doi: 10.1063/1.450106
*/

/*
 overlap integral
 unnormalized Cartesian Gaussian functions
 (a||b) = \int{dr \psi(r;\zeta_a,a,A) \psi(r;zeta_b,b,B)}
 %
 after integral of (s||s) = (\pi/\zeta)^(3/2) exp{-\xi(A-B)^2}
 where \zeta = \zeta_a + \zeta_b & \xi = \zeta_a * \zeta_b / \zeta
 and A B is the nuclear pos. so, we use R = |A-B|.
*/

/*
 un-normalized <s||s> integral
 $\zeta$ is the orbital exponent
 input: exponent of A,B. and the distance of A,B.
*/
 double 
 overlap_ssu(double dist, double za, double zb)
 {
    double zeta, xi;
    double k, ek, s;
    
    zeta = za + zb;
    xi = za * zb / zeta;
    k = pow(bM_PI/zeta, 3.0/2.0);
    ek = exp(-xi * dist * dist);
    s = k * ek;
    
    return s;
 }
 
 /*
 factorial of a number as n! or n!!
 itype == 1: n!
 itype == 2: n!!
*/
double
factorial (int n, int itype)
{
    int i;
    double s = 1.0;
    if (itype == 1 && n > 0) 
        for (i = 1; i <= n; i++)
            s *= i;
    else if (itype == 2 && n > 2)
        for(i = 1; i <= n; i += 2)
            s *= i;
    else
        s = 1.0;
    return s;
}


/*
    combination function use the formula : n!/(n-m)!/m!  
*/
double 
factorial_div(int n, int m)
{
    double s;

    s = factorial(n, 1) / factorial(m, 1) / factorial(n-m, 1);

    return s;
 }

int
total_angular_momentum(int ndx[3])
{
    int n;
    n = ndx[0] + ndx[1] + ndx[2];
    // printf("%12.6lf%12.6lf%12.6lf\n", ndx[0], ndx[1], ndx[2]);
    return n;    
}

/*
  caculate the normalization constant of one GTO, 
  argumented by n[3] index, and zeta exponent.
  Eq. 3 in JCP, 84, 3963 (1986)
  
  normal (\zeta,n) = {frac{2\zeta}{\pi}}^{3/4} (4\zeta)^{(n_x+n_y+n_z)/2} \times 
                     [(2n_x-1)!! (2n_y-1)!! (2n_z-1)!!]
  check: p_x orbital with \zeta = 0.1, we get 0.126739
*/
double 
norm_const_gto(int *ndx, double zeta)
{
    int nt;
    double k1, k2, k3, s;
 printf("%5d%5d%5d\n", ndx[0], ndx[1], ndx[2]);
    nt = total_angular_momentum(ndx);
    
    k1 = pow(2.0 * zeta / bM_PI, 0.75);
    k2 = pow(4.0 * zeta, nt * 0.5);
    k3 = factorial(2*ndx[0]-1, 2) * factorial(2*ndx[1]-1, 2) * factorial(2*ndx[2]-1, 2);
    printf("%5d%12.6lf%12.6lf%12.6lf\n", nt, k1, k2, k3);
    s = k1 * k2 * k3;
    
    return s;    
}


double 
norm_const_gaussian(int n[3],double zeta)
{
  double pre_c, up_c, fac_m, s1, s;

  pre_c = pow(2*zeta/bM_PI, 3.0e0/4.0);
  up_c = pow(4.0*zeta, (n[0]+n[1]+n[2]));
  fac_m = factorial(2*n[0]-1, 2) * factorial(2*n[1]-1, 2) * factorial(2*n[2]-1, 2);
  s1 = sqrt(up_c / fac_m);
  s = pre_c * s1;

  return s;
}


#undef PI

int main()
{
    double R = 1.0;
    double za = 0.1, zb = 0.1;
    double x;
    int i;
    double ndx[3];
    ndx[0] = 1; ndx[1] = 0; ndx[2] = 0;
     printf("%12.6lf%12.6lf%12.6lf\n", ndx[0], ndx[1], ndx[2]);

    // x = factorial(12, 1);
    x = norm_const_gto(ndx, za);

    printf("%12.6lf\n", x);
    // for (i = 0; i < 100; i++) {
        // x = overlap_ssu(R+0.1*i, za, zb);
        // printf("%12.8g\n", x);
    // }
    getchar();
    return 0;
}

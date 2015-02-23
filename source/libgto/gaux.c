
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
 int 
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




/*
 factorial of a number as n!
*/
double 
fac(int n)
{
    int i;
    double s = 1.0;
    if(n > 0)
        for(i = 1; i <= n; i++)
            s *= i;
    return s;
}


/*
 factorial of a number such as : n!!
*/
double 
fac_odd(int n)
{
    int i;
    double f = 1.0;

    if(n > 2) 
        for(i = 1; i <= n; i += 2)
            f *= i;
    return f;
}


 /*combination function use the formula : n!/(n-m)!/m!  */
double 
combine_fac(int n, int m)
{
    double s;

    s = fac(n) / fac(m) / fac(n-m);

    return s;
 }




/*
    caculation of Ni
    where Ni (n), standing for ith is meant to take the value of
    the i component of the angular momentum index n
    input: index & i & flag
    output: Ni - flag
*/
int 
angular_momentum_index(int index[3], int i, int flag)
{
   int x;
   
   x = index[i] - flag;
  
   return x;
}

/*
  this can determine the subshell number
  that is sum over index[3]
  */
  
int 
sum_index(int index[3])
  {
    int i;
    int n = 3, ndim = 3;
    for(i = 0; i < ndim; i++)
        n += index[i];

    return n;
  }
 
 

/*
caculate the 3D distance of coordinate of point A,B,and return the distance of the coordinate
*/

double 
point3d_dist2(double *ri,double *rj)
{
    int i;
    double x, y, z;
	double r2;
	// call vector norm routine
	r2 = 0.0;
	for(i = 0; i < n; i++) {
	    r2 += (rj[i] - ri[i]) * (rj[i] - ri[i]);
	}
	dist = sqrt(r2);
	
    return r2;
}

 /* 
   EXPLAIN:
   1. caculate square distance of A,B and the 3D distance is stored in the array *pC.
   2. double *zeta stored zetaA and zetaB 
   C = (\zeta_A * A + \zeta_B * B) / (\zeta_A + \zeta_B)
 */
 void
 weighted_site2(double *pA, double *pB, double *pC, double zetaA, double zetaB)
 {
    int i;
    int n3d = 3;
    for(i = 0; i < n3d; i++)
        pC[i] = (zetaA * pA[i] + zetaB * pB[i]) / (zetaA + zetaB);
    return;
 }

/*
   P = (\zeta_A * A + \zeta_B * B) / (\zeta_A + \zeta_B)
   Q = (\zeta_C * C + \zeta_D * D) / (\zeta_C + \zeta_D)
   W = (\zeta * P + \eta * Q) / (\zeta + \eta)
*/
void 
weighted_site4(double *pA, double *pB, double *pC, double *pD, double *pW,
               double zetaA, double zetaB, double zetaC, double zetaD)
{
    double z[2];
    double pointP[3], pointQ[3];
    double zeta, eta;
    
    zeta = zetaA + zetaB;
    eta = zetaC + zetaD;
    
    z[0] = zetaA; z[1] = zetaB;
    weighted_site2(pA, pB, pointP, z);
    
    z[0] = zetaC; z[1] = zetaD;
    weighted_site2(pC, pD, pointQ, z);
    
    z[0] = zeta; z[1] = eta;
    weighted_site2(pointP, pointQ, pW, z);
    
    return;
}


  /*
  1. INPUT index[3],and exponent, here it refers to GTO
  2.caculate AND return the normalization constant of GTO
  eq. 3; J. Chem. Phys. 84, 3963 (1986); doi: 10.1063/1.450106
  
  normal (\zeta,n) = {frac{2\zeta}{\pi}}^{3/4} (4\zeta)^{(n_x+n_y+n_z)/2} \times 
                     [(2n_x-1)!! (2n_y-1)!! (2n_z-1)!!]
  */
double 
normalize_gaussian(int n[3],double zeta)
{
  double pre_c, up_c, fac_m, s1, s;

  pre_c = pow(2*zeta/PI, 3.0e0/4.0);
  up_c = pow(4.0*zeta, n[0]+n[1]+n[2]);
  fac_m = fac_odd(2*n[0]-1) * fac_odd(2*n[1]-1) * fac_odd(2*n[2]-1);
  s1 = sqrt(up_c / fac_m);
  s = pre_c * s1;

  return s;
}


/*
    F_m(T) = \int_0^1{dt t^{2m} exp(-T t^2)}
    T = \rho (P-Q)^2
    the series formula...
    F_{m}(t) = e^(-t) \Sum{i=0}{\infinite}\frac{(2m-1)!!(2t)^i}{(2m+2i+1)!!}
*/
double 
Fm(int m, double t, double eps)
{
    int i;
    double s, v;
    
    i = 0;
    v = 0.0;
    s = 1.0 / (2.0*m+1.0);
    while(s > eps) {
        s = fac_odd(2*m-1) * pow(2.0*t,(float)i) / fac_odd(2*m+2*i+1);
        v += s;
        i++;
    }        
    v *= exp(-t);
    return v;
}















/*
  this struct is used for overlap, kinetic, eri's integral.
*/
typedef struct{

    double** gtovalue;
    double stovalue;

}Integral_t;


typedef struct{


}Cmath;


// atomic orbital integral between two atoms
typedef struct{
    int n_int;
    
}AtomicIntegral_t;


#include "vector.h"





#include "vector.h"

/*
  Reference:
  1.
  Efficient recursive computation of molecular integrals over Cartesian
  Gaussian functions
  S. Obara and A. Saika 
  J. Chem. Phys. 84, 3963 (1986); doi: 10.1063/1.450106
  2.
  General recurrence formulas for molecular integrals over Cartesian
  Gaussian functions
  S. Obara and A. Saika 
  J. Chem. Phys. 89, 1540 (1988); doi: 10.1063/1.455717 
  3.
  A method for twoelectron Gaussian integral and integral derivative
  evaluation using recurrence relations
  Martin HeadGordon and John A. Pople 
  J. Chem. Phys. 89, 5777 (1988); doi: 10.1063/1.455553 
*/


/*
  see ref. J. Chem. Phys. 84, 3963 (1986); 
  overlap for s-s type integral.
  Cartesian Gaussian functions
  formula:
  (s||s) = \left(\frac{\pi}{\zeta}\right) exp\{-\xi (\bf A - \bf B)^2\}
*/


void
set_integral_value(Integral_t *integral, int i, int j, int m, int n, double value)
{
    
}

double
get_integral_value()
{

}




static void 
ovlpu_ss(Integral_t **rs, GTOBasis_t *t, int i, int j, int m, int n)
{
    double expona, exponb, r2, r, s;  

    expona = t[i].expon * t[i].expon * t[i].expon_gto[m];
    exponb = t[j].expon * t[j].expon * t[j].expon_gto[n];

    // r2 = D3_AB2(t[i].coord, t[j].coord,ab);
    // r = sqrt(r2);
    r = points_distance(t[i].coord, t[j].coord, 3);
 
    s = ovlpu(expona, exponb, r);

    // rs[i][j].gtovalue[m][n] = s;

    // printf("\nrs[%d][%d].gtovalue[%d][%d] = %lf\n",i, j, m, n, rs[i][j].gtovalue[m][n]);
    return s;
 }



 /*
 <a+1i||s> 
 */
void ovlp_as(INTEGRALa **rs, STOa *t, int i, int j, int m, int n)
{
    int ni, nc2, pre_a1, pre_a2;

    double pi ,s1, s2, s, zeta1, zeta2, zeta, z;

    s = s1 = s2 = 0.0;
    z = 0.0;
    nc2=0;
    zeta1 = t[i].expon * t[i].expon * t[i].expon_gto[m];
    zeta2 = t[j].expon * t[j].expon * t[j].expon_gto[n];

    zeta = zeta1 + zeta2;


    ni = t[i].pre_d;

    pi = create_px_i(t[i].coord,t[j].coord,zeta1,zeta2,ni,'a');

    pre_a1=t[i].pre_a1;
    pre_a2=t[i].pre_a2;
    nc2=angular_momentum_index(t[pre_a1].index,ni,0);
 
    if(pi!=0.0) {
        s1=rs_seek(rs,pre_a1,j,m,n);
    }

    if(nc2>0) {
        z=1.0/2.0/zeta;
        s2=rs_seek(rs,pre_a2,j,m,n);
        s2=z*s2;
    }

    s=pi*s1+nc2*s2;
    rs[i][j].gtovalue[m][n]=s;
/*
 printf("%lfdddd",s);
 printf("\nrs[%d][%d].gtovalue[%d][%d]=%lf\n",i,j,m,n,rs[i][j].gtovalue[m][n]);
 */
}

/*
===============================<a||b+1i>==========================
*/


void ovlpu_ab(INTEGRALa **rs,STOa *t,int i,int j,int m,int n)
{
int ni,nc2,nc3,tmpi,tmpj;

double dc1,s1,s2,s3,s,zeta1,zeta2,zeta,z;


s1=s2=s3=s=0.0;
z=0.0;


  zeta1=t[i].expon*t[i].expon*t[i].expon_gto[m];
  zeta2=t[j].expon*t[j].expon*t[j].expon_gto[n];
  zeta=zeta1+zeta2;
  z=1.0/2.0/zeta;
  ni=t[j].pre_d;
  dc1=create_px_i(t[i].coord,t[j].coord,zeta1,zeta2,ni,'b');
  
  nc2=angular_momentum_index(t[j].index,ni,2);

  nc3=angular_momentum_index(t[i].index,ni,0);


  if(dc1!=0.0)
  {
   tmpj=t[j].pre_a1;

   s1=rs_seek(rs,i,tmpj,m,n);

  }
  if(nc2>0)
  {
  tmpj=t[j].pre_a2;

  s2=rs_seek(rs,i,tmpj,m,n);

  s2=s2*z; 


  }

  if(nc3>0)
  {
   tmpi=t[i].pre_a1;
   tmpj=t[j].pre_a1;

   s3=rs_seek(rs,tmpi,tmpj,m,n);

   s3=s3*z;

  }


  s=dc1*s1+nc2*s2+nc3*s3;

  rs[i][j].gtovalue[m][n]=s;

}


#undef PI

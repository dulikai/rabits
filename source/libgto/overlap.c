#include<stdio.h>
#include<math.h>

#include "int_tools.c"


/*
    Two-center overlap Integrals
    overlap integral between two GTOs
    \usepackage{bm}
    \usepackage{mathrsfs}
    \psi(\bm r;\zeta,\bm n,\bm R) = (x-R_x)^{n_x} (y-R_y)^{n_y} (z-R_z)^{n_z}
                        exp[-\zeta (r-R)^2]
    normalization:
    \mathcal {N} = {\frac{2\xi}{\pi}}^{3/4} (4\xi)^{(n_x+n_y+n_z)/2} \
                    [(2n_x-1)!!(2n_y-1)!!(2n_z-1)!!]^{-1/2}
*/





 void ovlpu_ss(INTEGRALa **rs,STOa *t,int i,int j,int m,int n)
 {
  double expona,exponb,r2,r,s,ab[3];
  

  expona=t[i].expon*t[i].expon*t[i].expon_gto[m];
  exponb=t[j].expon*t[j].expon*t[j].expon_gto[n];

  r2=D3_AB2(t[i].coord,t[j].coord,ab);

  r=sqrt(r2);
 
  s=ovlpu(expona,exponb,r);

  rs[i][j].gtovalue[m][n]=s;

 /*
  printf("\nrs[%d][%d].gtovalue[%d][%d]=%lf\n",i,j,m,n,rs[i][j].gtovalue[m][n]);
 */
 }

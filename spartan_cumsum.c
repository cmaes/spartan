#include "spartan.h"
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
double spartan_cumsum (spartan_index *p, spartan_index *c, spartan_index n)
{
  spartan_index i, nz = 0 ;
  double nz2 = 0 ;
  if (!p || !c) return (-1) ;     /* check inputs */
  for (i = 0 ; i < n ; i++)
    {
      p [i] = nz ;
      nz += c [i] ;
      nz2 += c [i] ;              /* also in double to avoid int overflow */
      c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
  p [n] = nz ;
  return (nz2) ;                  /* return sum (c [0..n-1]) */
}

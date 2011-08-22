#include "spartan.h"

/* compute d where d(j) = norm(A(:,j)) */
int spartan_norms(spartan_sparse *A, double *d, double dmin, double dmax){
  double alpha, aij, *x = A->x;
  spartan_index  p, j, *Ap = A->p, *Ai = A->i, n = A->n;
  if (!A || !d) return(-1);
  for (j = 0; j < n; j++){
    alpha = 0;
    for (p = Ap[j]; p < Ap[j+1]; p++){
      aij = x[Ai[p]];
      alpha += aij*aij;
    }
    alpha = sqrt(alpha);
    if (alpha < dmin) alpha = dmin; 
    if (alpha > dmax) alpha = dmax; 
    d[j] = alpha; 
  }
  return 0; 
}

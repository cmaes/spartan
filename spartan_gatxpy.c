#include "spartan.h" 

/* y = A'*x + y */ 
/* A: m x n => A': n x m => x : m x 1  => y = n x 1 */
int spartan_gatxpy (const spartan_sparse *A, const double *x, double *y){
  spartan_index p, j, n, *Ap, *Ai; 
  double *Ax; 
  if (!A || !x || !y) return(0);
  n = A->n; Ap = A->p; Ai = A->i;  Ax = A->x;
  for(j = 0; j < n; j++){
    for(p = Ap[j]; p < Ap[j+1]; p++){
      y[j] += Ax[p] * x[Ai[p]];
    }
  }
  return(1);
}

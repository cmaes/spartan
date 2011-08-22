#include <stdio.h>
#include "spartan.h"

int broydentridiagonalfunction(double *f, const double *x, spartan_index n){
  spartan_index i;
  double h = 0.5;
  
  f[0] = (3 - h*x[0])*x[0] - 2*x[1] + 1;
  for(i=1; i<=n-2; i++){
    f[i] = (3 - h*x[i])*x[i] - x[i-1] - 2*x[i+1] + 1;
  }
  f[n-1] = (3 - h*x[n-1])*x[n-1] - x[n-1] + 1;

  return(0);
}

int broydentridiagonaljacobian(spartan_sparse **J, const double *x,\
                               spartan_index n){
  spartan_index nnz, i;
  spartan_sparse *T;
  double h = 0.5;

  /* The number of nonzeros in the tridiagonal jacobian */
  nnz = 2 + 3*(n-2) + 2;  
  
  /* Allocate the Jacobian */ 
  T = spartan_spalloc(n,n,nnz,1,1);
  /* Check for out of memory error */ 
  if (T == NULL)
    return(-1);

  for (i = 0; i < n; i++){
    if (i > 0) spartan_entry(T,i,i-1,-1);
    
    spartan_entry(T, i, i, 3 - 2*h*x[i]);
    
    if (i < n-1) spartan_entry(T,i, i+1, -2);
  }

  /* convert sparse matrix from triplet format to
     compressed sparse column format */ 
  *J = spartan_compress(T);
  
  /* free the triplet form storage */
  T = spartan_spfree(T);

  return(0);
}

void printstring(const char *s){
  printf("%s\n",s);
}

int main(int argc, char **argv){

  spartan_index j, n = 1000;
  double *f, *x, opts[SPARTAN_OPTIONS]; 
  int status, i;
  /* Allocate space for f and x */ 

  f = malloc(sizeof(double)*n);
  x = malloc(sizeof(double)*n);

  if (!x || !f) return(-1);

  /* Setup default options */
  for(i=0; i<SPARTAN_OPTIONS; i++) opts[i] =0;
  /* Tell spartan to free the Jacobian after use */
  opts[SPARTAN_MEMMGMT] = 1;

  /* Set the print function for spartan output */
  spartan_setprintfunction(printstring);

  /* Set x = 0 */
  for (j=0;j<n;j++) x[j] = -1.;

  printf("running broyden tridiagonal demo n=%d\n",(int) n);
  status = spartan_solve(broydentridiagonalfunction, broydentridiagonaljacobian,\
                n, f, x, opts);
  printf("spartan_solve:%d\n",status);
  free(x); free(f);
  return(0);
}


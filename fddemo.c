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
  double h = 0.5;
  spartan_sparse *Jt;

  /* The number of nonzeros in the tridiagonal jacobian */
  nnz = 2 + 3*(n-2) + 2;  
  
  /* Allocate the Jacobian */ 
  Jt = spartan_spalloc(n,n,nnz,1,1);
  /* Check for out of memory error */ 
  if (Jt == NULL)
    return(-1);

  for (i = 0; i < n; i++){
    if (i > 0) spartan_entry(Jt,i,i-1,-1);
    
    spartan_entry(Jt, i, i, 3 - 2*h*x[i]);
    
    if (i < n-1) spartan_entry(Jt,i, i+1, -2);
  }

  /* convert sparse matrix from triplet format to
     compressed sparse column format */ 
  *J = spartan_compress(Jt);
  
  spartan_spfree(Jt);
  return(0);
}

void printstring(const char *s){
  printf("%s\n",s);
}

int main(int argc, char **argv){

  spartan_index j, k, n = 1024, nnz, *coloring, numcolors, mincolors;
  double *f, *x, opts[SPARTAN_OPTIONS]; 
  spartan_sparse *Jpat_triplet, *Jpat_compress, *J_fd, *J_usr;
  spartan_jacobian sparsefdjacobian;
  int status, i;
  /* Allocate space for f and x */ 

  f = malloc(sizeof(double)*n);
  x = malloc(sizeof(double)*n);

  if (!x || !f) return(-1);

  /* Setup default options */
  for(i=0; i<SPARTAN_OPTIONS; i++) opts[i] =0;
  
  /* Setup the print function */
  spartan_setprintfunction(printstring);

  /* Set x = -e */
  for (j=0;j<n;j++) x[j] = -1.;

  /* Construct the sparsity pattern of the Jacobian */ 
  /* The number of nonzeros in the tridiagonal jacobian */
  nnz = 2 + 3*(n-2) + 2;  
    
  /* Allocate the Jacobian pattern */ 
  Jpat_triplet = spartan_spalloc(n,n,nnz,1,1);
  /* Check for out of memory error */ 
  if (Jpat_triplet == NULL)
    return(-1);
    
  /* Fill in the sparsity pattern of the Jacobian */
  for (i = 0; i < n; i++){
    if (i > 0) spartan_entry(Jpat_triplet,i,i-1,1);
      
    spartan_entry(Jpat_triplet, i, i, 1);
    
    if (i < n-1) spartan_entry(Jpat_triplet,i, i+1, 1);
  }
  
  /* Allocate the column coloring */
  coloring = spartan_malloc(n,sizeof(spartan_index));

  /* Compute a column coloring of the Jacobian pattern */ 
  spartan_columncolor(coloring,&numcolors,&mincolors,Jpat_triplet);

  printf("computed a coloring with %d colors. Lower bound:%d Max:%d\n", 
         (int) numcolors, (int) mincolors, (int) n);

  /* Convert the Jacobian pattern to compressed-sparse column format */ 
  Jpat_compress = spartan_compress(Jpat_triplet);

  /* free triplet format */
  spartan_spfree(Jpat_triplet);

  /* Get the function to compute the sparse Jacobian via finite differences */
  sparsefdjacobian = spartan_initsfd(broydentridiagonalfunction, Jpat_compress, coloring, numcolors);

  printf("running broyden tridiagonal (sparse finite differences) demo n=%d\n",(int) n);
  status = spartan_solve(broydentridiagonalfunction, sparsefdjacobian,\
                n, f, x, opts);
  printf("spartan_solve:%d\n",status);


  spartan_spfree(Jpat_compress);
  spartan_free(coloring);
  free(x);
  free(f);

  return(0);
}


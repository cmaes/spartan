#include "spartan.h" 

/* [coloring,numcolors,mincolors] = columncolor(J)

   Input:   J           A sparse m x n Jacobian matrix in triplet form. 

   Output:  coloring    A coloring of the columns of J. Of length n. 
            numcolors   The number of colors used in the coloring 
            mincolors   A lower bound on the minimum number of colors
                        required to color the columns of Jacobian

   columncolor calls DSM from TOMS Algorithm 618 to compute 
   a column coloring of a sparse Jacobian */ 

/* DSM won't work with 64-bit indicies so we only define the routine 
   if LP64 is not definied. This will lead to a link error for users
   on 64 bit platforms. */
#ifndef LP64 
int spartan_columncolor(int *coloring, int *numcolors,
                int *mincolors, const spartan_sparse *J){

  /* A work array */
  int *w, k, *Jp, *Ji;

  /* The following variables are used for the call to DSM */
  int m, n, nnz;
  int *indrow, *indcol;
  int info; 
  int *ipntr;
  int *jpntr;
  int *iworkarray; 
  int lenwrk;

  /* Get m and n and column and row pointers*/ 
  m = J->m; n = J->n; Jp = J->p; Ji = J->i;

  /* Get the number of nonzeros in J */
  nnz = J->nzmax; 

  /* Calculate the size of DSMs workspace */ 
  lenwrk = m; 
  if (lenwrk < 6*n) 
    lenwrk = 6*n;

  /* Allocate memory */ 
  w = spartan_malloc(2*nnz + (m+1) + (n+1) + lenwrk, sizeof(int));

  if (!w) return (-1);

  /* Get pointers into the work array */
  indrow     = w;
  indcol     = w + nnz;
  ipntr      = w + nnz + nnz; 
  jpntr      = w + nnz + nnz + (m+1);
  iworkarray = w + nnz + nnz + (m+1) + (n+1);

  /* Convert C indices to Fortan indicies */
  for(k=0; k<nnz; k++){
    indcol[k] = Jp[k]+1;
    indrow[k] = Ji[k]+1;
  }

  /* Call DSM to perform the actual work */
  dsm_(&m,&n,&nnz,indrow,indcol,coloring,numcolors,mincolors,&info,ipntr,jpntr,iworkarray,&lenwrk);

  /* Convert Fortan indices to C indicies */
  for (k=0; k<n; k++) coloring[k]--;

  w = spartan_free(w);

  return (info);
}
#endif

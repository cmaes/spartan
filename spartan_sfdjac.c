#include "spartan.h"
/*
  J = spartan_sfdjac(userfun,x,J,coloring) Computes a finite-difference 
  approximation of a sparse Jacobian 
 
  Input:  userfun      A user-defined function of the form 
                       f = userfun(x) 
                       which returns f(x), a vector of length m
 
          x            A vector of length n, at which to evaluate J(x) 
 
          coloring     A coloring of the columns of the Jacobian such 
                       that coloring(j) = k implies that the jth column 
                       of J is colored k.  
                        
                       A coloring is a partition of the columns of 
                       the Jacobian into distinct groups (or colors) such 
                       that no column of the same color has a nonzero entry 
                       in the same row.
  
                       All columns of the same color can be estimated via 
                       a single finite-difference. Thus, the number of
                       of function evaluations required is equal to the number
                       of colors plus one.

          numcolors    The number of colors in the coloring

  Input/Output:
          J           
                       A m x n matrix that contains the sparsity pattern 
                       of the Jacobian. J should be compressed-column format, 
                       storage for the nonzero values should be allocated. 
                       The nonzero values of J are altered.
*/

spartan_sparse *spartan_sfdjac
( spartan_function usrfun,
  const double *x, 
  spartan_sparse *J, 
  spartan_index *coloring, 
  spartan_index numcolors
)
{
  spartan_index n, j, k, color; 
  spartan_index *Ji, *Jp;
  double *w, *f, *xp, *fp;
  double *Jx;

  /* For now fix epsilon as the sqrt of 
     machine precion. The user should be 
     allowed to pass in a vector of epsilons */
  double epsilon = sqrt(DBL_EPSILON);
  int err;

  n = J->n; Ji = J->i; Jp = J->p; Jx = J->x;
  
  /* allocate work arrays */
  w  = spartan_malloc(3*n,sizeof(double));
  if(!w) return (NULL);
  f  = w;
  xp = w+n;
  fp = w+2*n;


  /* compute f = F(x) */ 
  err = usrfun(f,x,n);
  if (err){ spartan_free(w);  return (NULL); };
  
  for (color=0; color<numcolors; color++){
    /* Construct xp */ 
    /* xp[j] = x[j] + epsilon if column j is colored color */ 
    /* otherwise xp[j] = x[j] */ 
    for(j=0; j<n; j++){
      if (coloring[j] == color)
        xp[j] = x[j] + epsilon;
      else 
        xp[j] = x[j]; 
    }

    /* Compute fp = userfun(xp) */ 
    err = usrfun(fp,xp,n);
    if (err){ spartan_free(w);  return (NULL); };
    
    /* fp ~= f + J*(epsilon*p) */

    /* fp <- (fp - f)/epsilon */ 
    for(j=0; j<n; j++)
      fp[j] = (fp[j] - f[j])/epsilon;
    
    /* place the nonzeros in fp in J */
    for(j=0; j<n; j++){
      if (coloring[j] == color){
        for(k = Jp[j]; k< Jp[j+1]; k++) /* Loop over nonzeros in col j */ 
          Jx[k] = fp[Ji[k]]; /* Set J(i,j) = fp(i) */   
      }
    }    
  }

  return (spartan_done(J,w,NULL,1));
}




/*Global variables to store the nonzero pattern of the Jacobian
  and a column coloring */
spartan_sparse *spartan_Jpattern;  
spartan_index  *spartan_colcoloring;
spartan_index  spartan_numcolors;
spartan_function spartan_fdusr;

/* prototype for return to compute Jacobian */
int spartan_fdjacobian(spartan_sparse **J, const double *x, spartan_index n);

spartan_jacobian spartan_initsfd(spartan_function usrfun, spartan_sparse *Jpattern,
                                 spartan_index *colcoloring, spartan_index numcolors){
  spartan_fdusr       = usrfun;
  spartan_Jpattern    = Jpattern;
  spartan_colcoloring = colcoloring;
  spartan_numcolors   = numcolors;
  return (spartan_fdjacobian);
}

int spartan_fdjacobian(spartan_sparse **J, const double *x, spartan_index n){
  spartan_index j,k;
  *J = spartan_sfdjac(spartan_fdusr,x,spartan_Jpattern,spartan_colcoloring, spartan_numcolors);
  if (!*J) return(-1);
  return (0);
}

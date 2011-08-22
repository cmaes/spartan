#include "spartan.h"
#include <string.h> /* for memset */

/* alpha = ||x|| */ 
double spartan_dnrm2(const double *x, const spartan_index n){
  return (sqrt(spartan_ddot(x,x,n)));
}

/* alpha = x'*y */ 
double spartan_ddot (const double *x, const double *y,  const spartan_index n){
  double alpha = 0; 
  spartan_index j;
  for (j = 0; j < n; j++)
    alpha += x[j]*y[j];
  return(alpha);
}
/* x <- alpha*x */ 
int    spartan_dscal(const double alpha, double *x, const spartan_index n){
  spartan_index j;
  if (!x) return(-1);
  for (j = 0; j < n; j++)
    x[j] *= alpha; 
  return(0);
}

/* y <- alpha*x + y */  
int    spartan_daxpy(const double alpha, const double *x, double *y, const spartan_index n){
  spartan_index j;
  if (!x || !y) return(0);
  for(j=0;j<n;j++)
    y[j] += alpha*x[j];
  return(0);
}

/* y <- x */ 
int    spartan_dcopy(const double *x, double *y, const spartan_index n){
  if (!x || !y) return(-1);
  memcpy(y,x,sizeof(double)*n);
  return(0);
}

/* Non-standard (but useful) BLAS */

/*  z <- x.*y */ 
int    spartan_dmult(const double *x, const double *y, double *z, 
                     const spartan_index n){
  if (!x || !y || !z) return(-1);
  spartan_index j;
  for (j=0; j < n; j++)
    z[j] = x[j]*y[j]; 
  return 0;
}

/* || f ||_inf */
double spartan_norminf(const double *f, spartan_index n){
  spartan_index j;
  double max = 0;
  double absfj; 
  if (!f) return (0);
  for(j=0; j<n; j++){
    absfj = fabs(f[j]);
    if (absfj > max) max = absfj;
  }
  return max;
}

/* x <- 0 */
int spartan_zero(double *x, spartan_index n){
  spartan_index j;
  if (!x) return(-1);
  for(j=0; j<n; j++) x[j] = 0.;
  return(0);
}

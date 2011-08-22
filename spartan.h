#ifndef _SPARTAN_H
#define _SPARTAN_H
#include <stdlib.h> 
#include <math.h> 
#include <float.h>
#include <string.h>

#ifdef MATLAB_MEX 
#ifndef MX_COMPAT_32
#ifndef LP64
#define LP64
#endif 
#endif 
#endif 

/* This typedef allows us to deal with largeArrayDims */
#ifdef LP64
typedef long long  spartan_index;
#else 
typedef int        spartan_index;
#endif 

typedef struct spartan_sp {
  spartan_index nzmax; /* maximum number of entries */
  spartan_index  m;    /* number of rows */ 
  spartan_index  n;    /* number of cols */ 
  spartan_index *p;    /* column pointers (size n+1) or col indices (size nzmax)*/ 
  spartan_index *i;    /* row indices, size nzmax */
  double        *x;    /* numerical values, size nzmax */  
  spartan_index nz;    /* # of entries in triplet form, -1 for compressed-col*/ 
} spartan_sparse;

#define SP_CSC(A) (A && (A->nz == -1))
#define SP_TRIPLET(A) (A && (A->nz >= 0))

/* sparse utility routines */ 
spartan_sparse *spartan_spalloc (spartan_index m, spartan_index n,\
                                 spartan_index nzmax,int values, int triplet);
int spartan_sprealloc(spartan_sparse *A, spartan_index nzmax);
spartan_sparse *spartan_spfree(spartan_sparse *A);
spartan_sparse *spartan_done (spartan_sparse *C, void *w, void *x, int ok);
int spartan_entry(spartan_sparse *A, spartan_index i, spartan_index j, double x);
spartan_sparse *spartan_compress (const spartan_sparse *T);
double spartan_cumsum(spartan_index *p, spartan_index *c, spartan_index n);
int spartan_norms(spartan_sparse *J, double *d, double dmin, double dmax);

/* sparse matrix-vector products */ 
int spartan_gaxpy (const spartan_sparse *A, const double *x, double *y);
int spartan_gatxpy(const spartan_sparse *A, const double *x, double *y);

/* vector operations (Level-1 BLAS)*/ 
double spartan_dnrm2(const double *x, const spartan_index n);
double spartan_ddot (const double *x, const double *y,  const spartan_index n);
int    spartan_dscal(const double a, double *x, const spartan_index n);
int    spartan_daxpy(const double a, const double *x, double *y, const spartan_index n);
int    spartan_dcopy(const double *x, double *y, const spartan_index n);

/* vector (non BLAS operations) */ 
int    spartan_dmult(const double *x, const double *y, double *z, 
                     const spartan_index n);
double spartan_norminf(const double *f, spartan_index n);
int    spartan_zero(double *x, spartan_index n);

/*  Linear solve */ 
#define SPARTAN_SINGULAR_JACOBIAN -1
int spartan_linsolve(const spartan_sparse *A, const double *b, double *x);

/* print output */ 
typedef void (*spartan_printfunc)(const char *string);
void spartan_print(const char *string);
void spartan_setprintfunction(spartan_printfunc printfunction);

/* Utilities */ 
void *spartan_calloc(spartan_index n, size_t size);
void *spartan_malloc(spartan_index n, size_t size);
void *spartan_realloc (void *p, spartan_index n, size_t size, int *ok);
void *spartan_free(void *ptr);


/* Computational routines */ 
int spartan_dogleg(const double *f, const spartan_sparse *J, const double Delta,
                   const double *d, double *p);

/* Function pointer types */
typedef int (*spartan_function)(double *f, const double *x, spartan_index n);
typedef int (*spartan_jacobian)(spartan_sparse **J, const double *x, spartan_index n);

/* Finite Difference init routine */
spartan_jacobian spartan_initsfd(spartan_function usrfun, spartan_sparse *Jpattern,
                                 spartan_index *colcoloring, spartan_index numcolors);


/* Main routine */ 
int spartan_solve(spartan_function usrfun,
                  spartan_jacobian usrjac,
                  spartan_index n, double *fval, double *x, 
                  double *opts);

/* Options */ 
#define SPARTAN_DELTAMAX 0 
#define SPARTAN_DELTA0   1 
#define SPARTAN_ETA      2
#define SPARTAN_TOL      3 
#define SPARTAN_JACSTOL  4
#define SPARTAN_ITERMAX  5
#define SPARTAN_DMIN     6
#define SPARTAN_DMAX     7 
#define SPARTAN_SCALING  8 
#define SPARTAN_MEMMGMT  9 

/* The number of options */
#define SPARTAN_OPTIONS  10      

/* epsilon for ieee double precision arithmetic */ 
#define MACHINE_EPS  DBL_EPSILON

#define SPARTAN_DELTAMAX_VALUE   1e2
#define SPARTAN_DELTA0_VALUE     1e0
#define SPARTAN_ETA_VALUE        1e-3
#define SPARTAN_TOL_VALUE        1e-6
#define SPARTAN_JACSTOL_VALUE    1e-12
#define SPARTAN_ITERMAX_VALUE    50
#define SPARTAN_DMIN_VALUE       MACHINE_EPS
#define SPARTAN_DMAX_VALUE       1/MACHINE_EPS
#define SPARTAN_SCALING_VALUE    0
#define SPARTAN_MEMMGMT_VALUE    0 


/* Error codes */ 
#define SPARTAN_OK               0 
#define SPARTAN_MEMORY_FAILURE  -1 
#define SPARTAN_USER_ERROR      -2
#define SPARTAN_LOCAL_MIN       -3
#define SPARTAN_ITER_LIMIT      -4

/* Macros */
#define SP_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SP_MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include "spartan.h"

/* allocate a sparse matrix (triplet form or compressed-column form) */
spartan_sparse *spartan_spalloc (spartan_index m, spartan_index n, spartan_index nzmax,
                                 int values, int triplet){
  spartan_sparse *A = spartan_calloc (1, sizeof (spartan_sparse)) ;    /* allocate the cs struct */
  if (!A) return (NULL) ;                 /* out of memory */
  A->m = m ;                              /* define dimensions and nzmax */
  A->n = n ;
  A->nzmax = nzmax = SP_MAX (nzmax, 1) ;
  A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
  A->p = spartan_malloc (triplet ? nzmax : n+1, sizeof (spartan_index)) ;
  A->i = spartan_malloc (nzmax, sizeof (spartan_index)) ;
  A->x = values ? spartan_malloc (nzmax, sizeof (double)) : NULL ;
  return ((!A->p || !A->i || (values && !A->x)) ? spartan_spfree (A) : A) ;
}


/* change the max # of entries sparse matrix */
int spartan_sprealloc (spartan_sparse *A, spartan_index nzmax)
{
  int ok, oki, okj = 1, okx = 1 ;
  if (!A) return (0) ;
  if (nzmax <= 0) nzmax = (SP_CSC (A)) ? (A->p [A->n]) : A->nz ;
  A->i = spartan_realloc (A->i, nzmax, sizeof (spartan_index), &oki) ;
  if (SP_TRIPLET (A)) A->p = spartan_realloc (A->p, nzmax, sizeof (spartan_index), &okj) ;
  if (A->x) A->x = spartan_realloc (A->x, nzmax, sizeof (double), &okx) ;
  ok = (oki && okj && okx) ;
  if (ok) A->nzmax = nzmax ;
  return (ok) ;
}

/* free a sparse matrix */
spartan_sparse *spartan_spfree (spartan_sparse *A)
{
  if (!A) return (NULL) ;     /* do nothing if A already NULL */
  spartan_free (A->p) ;
  spartan_free (A->i) ;
  spartan_free (A->x) ;
  return (spartan_free (A)) ;      /* free the cs struct and return NULL */
}


/* free workspace and return a sparse matrix result */
spartan_sparse *spartan_done (spartan_sparse *C, void *w, void *x, int ok)
{
  spartan_free (w) ;                       /* free workspace */
  spartan_free (x) ;
  return (ok ? C : spartan_spfree (C)) ;   /* return result if OK, else free it */
}

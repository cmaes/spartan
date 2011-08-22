#include "spartan.h"
/* C = compressed-column form of a triplet matrix T */
spartan_sparse *spartan_compress (const spartan_sparse *T)
{
  spartan_index m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
  double *Cx, *Tx ;
  spartan_sparse *C ;
  if (!SP_TRIPLET (T)) return (NULL) ;                 /* check inputs */
  m = T->m ; n = T->n ; Ti = T->i ; Tj = T->p ; Tx = T->x ; nz = T->nz ;
  C = spartan_spalloc (m, n, nz, Tx != NULL, 0) ;      /* allocate result */
  w = spartan_calloc (n, sizeof (spartan_index)) ;     /* get workspace */
  if (!C || !w) return (spartan_done (C, w, NULL, 0)) ;/* out of memory */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;            /* column counts */
  spartan_cumsum (Cp, w, n) ;                          /* column pointers */
  for (k = 0 ; k < nz ; k++)
    {
      Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
      if (Cx) Cx [p] = Tx [k] ;
    }
  return (spartan_done (C, w, NULL, 1)) ;      /* success; free w and return C */
}

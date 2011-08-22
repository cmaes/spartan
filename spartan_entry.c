#include "spartan.h"
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
int spartan_entry (spartan_sparse *T, spartan_index i, spartan_index j, double x)
{
  if (!SP_TRIPLET (T) || i < 0 || j < 0) return (0) ;     /* check inputs */
  if (T->nz >= T->nzmax && !spartan_sprealloc (T,2*(T->nzmax))) return (0) ;
  if (T->x) T->x [T->nz] = x ;
  T->i [T->nz] = i ;
  T->p [T->nz++] = j ;
  T->m = SP_MAX (T->m, i+1) ;
  T->n = SP_MAX (T->n, j+1) ;
  return (1) ;
}

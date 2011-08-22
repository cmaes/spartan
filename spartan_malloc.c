#include "spartan.h" 

/* There seems to be an issue with using the mx memory 
   management options. No errors are returned but Matlab
   hangs and will not execute another command 
#ifdef MATLAB_MEX_FILE
#include "mex.h" 
#define malloc mxMalloc
#define free   mxFree 
#define calloc mxCalloc
#endif 
*/

/* wrapper for malloc */
void *spartan_malloc (spartan_index n, size_t size){
  return (malloc (SP_MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *spartan_calloc (spartan_index n, size_t size){
  return (calloc (SP_MAX (n,1), size)) ;
}

/* wrapper for free */
void *spartan_free (void *p){
  if (p) free (p) ;       /* free p if it is not already NULL */
  return (NULL) ;         /* return NULL to simplify the use of spartan_free */
}

/* wrapper for realloc */
void *spartan_realloc (void *p, spartan_index n, size_t size, int *ok){
  void *pnew ;
  pnew = realloc (p, SP_MAX (n,1) * size) ; /* realloc the block */
  *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
  return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

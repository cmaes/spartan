#include "spartan.h" 

/* Perform the following test to get UMFPACKs header file 
   from the suitesparse directory inside of Scilab */
#ifndef SPARTAN_SCILAB  
#include "umfpack.h"
#else
#include <suitesparse/umfpack.h>
#endif 

#include <stdio.h> /* for sprintf */ 

/*  Ax = b */ 
#ifdef LP64
int spartan_linsolve(const spartan_sparse *A, const double *b, double *x){
  /* Using UMFPACKs sparse LU solve */ 
  UF_long status=UMFPACK_OK;
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric; 
  

  umfpack_dl_defaults(Control);
  status = umfpack_dl_symbolic(A->m,A->n,A->p,A->i,A->x,&Symbolic,Control,Info);
  if (status != UMFPACK_OK){
    char buf[1024]; 
    sprintf(buf,"symbolic failed:%ld", status); 
    umfpack_dl_report_status(Control,status);
    Control[UMFPACK_PRL] = 3;
    umfpack_dl_report_matrix(A->m,A->n,A->p,A->i,A->x,1,Control);
    spartan_print(buf);
    umfpack_dl_free_symbolic(&Symbolic);
    return (int) status; 
  } 

  status = umfpack_dl_numeric (A->p,A->i,A->x,Symbolic,&Numeric,Control,Info);
  if (status != UMFPACK_OK ){
    char buf[1024]; 
    UF_long stat2;
    umfpack_dl_report_status(Control,status);
    stat2 = umfpack_dl_report_matrix (A->m, A->n, A->p, A->i, A->x, 1, Control) ;
    sprintf(buf,"numeric failed:%ld %ld", status,stat2); 
    spartan_print(buf);
    umfpack_dl_free_symbolic(&Symbolic);
    umfpack_dl_free_numeric(&Numeric);
    return SPARTAN_SINGULAR_JACOBIAN; 
  } 

  umfpack_dl_free_symbolic(&Symbolic);
  status = umfpack_dl_solve(UMFPACK_A,A->p,A->i,A->x,x,b,Numeric,Control,Info);
    if (status != UMFPACK_OK){
    char buf[1024]; 
    sprintf(buf,"solve failed:%ld", status); 
    spartan_print(buf);
    umfpack_dl_free_symbolic(&Symbolic);
    umfpack_dl_free_numeric(&Numeric);
    return (int) status; 
  } 


  umfpack_dl_free_numeric(&Numeric);

  return(0);
}
#else 

int spartan_linsolve(const spartan_sparse *A, const double *b, double *x){
  /* Using UMFPACKs sparse LU solve */ 
  int status=UMFPACK_OK;
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric; 
  
  umfpack_di_defaults(Control);
  status = umfpack_di_symbolic(A->m,A->n,A->p,A->i,A->x,&Symbolic,Control,Info);
  if (status != UMFPACK_OK){
    char buf[1024]; 
    sprintf(buf,"symbolic failed:%d", status); 
    spartan_print(buf);
    umfpack_di_free_symbolic(&Symbolic);
    return (int) status; 
  } 

  status = umfpack_di_numeric (A->p,A->i,A->x,Symbolic,&Numeric,Control,Info);
  if (status != UMFPACK_OK ){
    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);
    return SPARTAN_SINGULAR_JACOBIAN; 
  } 

  umfpack_di_free_symbolic(&Symbolic);
  status = umfpack_di_solve(UMFPACK_A,A->p,A->i,A->x,x,b,Numeric,Control,Info);
    if (status != UMFPACK_OK){
    char buf[1024]; 
    sprintf(buf,"solve failed:%d", status); 
    spartan_print(buf);
    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);
    return (int) status; 
  } 


  umfpack_di_free_numeric(&Numeric);

  return(0);
}

#endif

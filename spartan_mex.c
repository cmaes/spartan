#include "spartan.h"
#include "mex.h" 
#include "umfpack.h"

int spartan_mex_check(int m, int n, int values, int sparse, int square,const mxArray *A){
  spartan_index mm = mxGetM(A), nn = mxGetN(A);
  if (values){
    if(mxIsComplex(A)) mexErrMsgTxt("matrix must be real");
    if(!mxIsDouble(A)) mexErrMsgTxt("matrix must be double");
  }
  if (sparse && !mxIsSparse(A)) mexErrMsgTxt("matrix must be sparse");
  if (!sparse && mxIsSparse(A)) mexErrMsgTxt("matrix must be full");
  if (m >= 0 && m != mm) mexErrMsgTxt("wrong dimension");
  if (n >= 0 && n != nn) mexErrMsgTxt("wrong dimension"); 
  if (square && mm != nn) mexErrMsgTxt("matrix must be square");
  return(0);
}

spartan_sparse *spartan_mex_get_sparse(spartan_sparse *A, int values, int square, 
                                      const mxArray *Amatlab){
  spartan_mex_check(-1,-1,values,1,square,Amatlab);
  A->m = mxGetM(Amatlab);
  A->n = mxGetN(Amatlab);
  A->p = mxGetJc(Amatlab);
  A->i = mxGetIr(Amatlab);
  A->x = values ? mxGetPr(Amatlab) : NULL; 
  A->nz = mxGetNzmax(Amatlab);
  return(A);
}

double *spartan_mex_get_double(spartan_index n, const mxArray *x){
  spartan_mex_check(n,1,1,0,0,x);
  return (mxGetPr(x));
}

void mex_print(const char *string){
  mexPrintf("%s\n",string); 
}

/* p = dogleg(f,J,Delta,d) */ 
void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]){
  
  spartan_index n; 
  spartan_sparse *J, Jmatrix; 
  double *f, *d, *Deltaptr, Delta, *p;
  const mxArray *frhs, *Jrhs, *Deltarhs, *drhs;
  mxArray *pl;

  if ( sizeof(spartan_index) != sizeof(UF_long) || 
       sizeof(spartan_index) != sizeof(mwIndex) )
    mexErrMsgTxt("Error incompatiables index sizes");
  
  /* Check number of input arguments */ 
  if (nrhs != 4)
    mexErrMsgTxt("spartan_mex(f,J,Delta,d)\n");
  
  frhs = prhs[0]; Jrhs = prhs[1]; Deltarhs = prhs[2]; drhs = prhs[3];
  pl   = plhs[0];
  
  n = mxGetM(frhs);

  /* Get input arguments */ 
  f = spartan_mex_get_double(n, frhs);
  spartan_mex_check(n, n, 1, 1, 1, Jrhs);
  J = spartan_mex_get_sparse(&Jmatrix, 1, 1, Jrhs);
  Deltaptr = spartan_mex_get_double(1, Deltarhs);
  Delta = Deltaptr[0];
  d = spartan_mex_get_double(n, drhs);

  /* Create output vector */ 
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  p  = spartan_mex_get_double(n,plhs[0]);

  /* For output */ 
  spartan_setprintfunction(mex_print);
  spartan_dogleg(f, J, Delta, d, p);
}

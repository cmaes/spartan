#include <stdio.h> 
#include "spartan.h"

double spartan_get_opt(double *opts, int option){
  double spartan_default_opts[]= { SPARTAN_DELTAMAX_VALUE,
                                   SPARTAN_DELTA0_VALUE,
                                   SPARTAN_ETA_VALUE,
                                   SPARTAN_TOL_VALUE,
                                   SPARTAN_JACSTOL_VALUE,
                                   SPARTAN_ITERMAX_VALUE,
                                   SPARTAN_DMIN_VALUE,
                                   SPARTAN_DMAX_VALUE,
                                   SPARTAN_SCALING_VALUE,
                                   SPARTAN_MEMMGMT_VALUE
                              };
  if (opts == NULL || opts[option] == 0.)
    return(spartan_default_opts[option]);
  else
    return(opts[option]);
}


int spartan_solve(int (*usrfun)(double *f, const double *x,spartan_index n),
                  int (*usrjac)(spartan_sparse **J, const double *x,spartan_index n),
                  spartan_index n, double *f, double *x, 
                  double *opts){
  
  /* options */
  double DeltaMax, Delta, eta, tol, jacstol, dmin, dmax; 
  spartan_index itermax;
  int scaling;
  int memmgmt;

  /* buffer for output */
  char output[1024];

  /* algorithm quantities */
  spartan_index k=0;
  int converged=0;
  double *work, *d, *xp, *fp, *fpred, *p, *grad, *w;
  double normDp, Delta_old, rho, obj, obj_true, obj_model, actual, pred, finf;
  spartan_index j; 
  spartan_sparse *J;
  int status; 

  DeltaMax = spartan_get_opt(opts,SPARTAN_DELTAMAX);
  Delta    = spartan_get_opt(opts,SPARTAN_DELTA0);
  eta      = spartan_get_opt(opts,SPARTAN_ETA);
  tol      = spartan_get_opt(opts,SPARTAN_TOL);
  jacstol  = spartan_get_opt(opts,SPARTAN_JACSTOL);
  dmin     = spartan_get_opt(opts,SPARTAN_DMIN);
  dmax     = spartan_get_opt(opts,SPARTAN_DMAX);

  scaling  = (int)   spartan_get_opt(opts,SPARTAN_SCALING);
  itermax  = (spartan_index) spartan_get_opt(opts,SPARTAN_ITERMAX);
  memmgmt  = (int)   spartan_get_opt(opts,SPARTAN_MEMMGMT);

  /* allocate workspace memory */
  work  = spartan_calloc(7*n,sizeof(double));
 
  if(!work) return SPARTAN_MEMORY_FAILURE;

  /* Create pointers into the workspace */
  d     = work;
  xp    = work+n;
  fp    = work+2*n;
  fpred = work+3*n;
  p     = work+4*n;
  grad  = work+5*n;
  w     = work+6*n;

  /* Call the users function */ 
  status = usrfun(f,x,n);
  if (status){
    spartan_print("error on first call");
    /* Clean up */ 
    spartan_free(work);
    return SPARTAN_USER_ERROR;
  }

  /* Compute the merit (objective) function */ 
  obj = spartan_ddot(f,f,n);
  
  /* Call the users Jacobian */
  status = usrjac(&J,x,n);
  if(status){
    spartan_print("error on first jac call");
    /* Clean up */ 
    spartan_free(work);
    return SPARTAN_USER_ERROR;
  }

  if(scaling) spartan_norms(J,d,dmin,dmax);
  else for(j=0; j<n; j++) d[j] = 1;
  
  /* Create the header */
  snprintf(output,sizeof(char)*1024,"SPARTAN: SPARse Trust-region Algorithm for Nonlinear equations");
  spartan_print(output);
  
  snprintf(output,sizeof(char)*1024,"itn   ||f||    rho     ||Dp||    Delta    err");
  spartan_print(output);

  finf = spartan_norminf(f,n);
  snprintf(output,sizeof(char)*1024,"%3d %8.1e                   %8.1e",\
           (int) k,finf,Delta);
  spartan_print(output);

  while( !converged && k < itermax){
    spartan_dogleg(f,J,Delta,d,p);

    spartan_dmult(d,p,w,n);  /* w <- d.*p */
    normDp = spartan_dnrm2(w,n);      /* || D*p || */
    Delta_old = Delta;
    
    /* Compute the function value predicted at this step */
    spartan_dcopy(f,fpred,n);   /* fpred <- f */
    spartan_gaxpy(J,p,fpred);   /* fpred = f + J*p */
    obj_model = spartan_ddot(fpred,fpred,n);
    
    /* Compute the function value at this step */
    spartan_dcopy(x,xp,n);      /* xp <- x */ 
    spartan_daxpy(1.,p,xp,n);   /* xp = p + x */
    status = usrfun(fp,xp,n);
    if (status){
      /* Clean up */ 
      spartan_free(work);
      return SPARTAN_USER_ERROR;
    }
    obj_true  = spartan_ddot(fp,fp,n);
    
    actual = obj - obj_true;
    pred   = obj - obj_model;

    sprintf(output,"obj:%e obj_true:%e obj_model:%e",obj,obj_true,obj_model);
    /*spartan_print(output);*/

    rho = actual/pred; 

    if(rho < 0.25)
      Delta = 0.25*Delta;
    else{
      if (rho > 0.75 && fabs(normDp - Delta) < 1e-8){
        Delta = SP_MIN(2*Delta,DeltaMax);
      }
    }

    if (rho > eta){
      spartan_daxpy(1,p,x,n);  /* x <- p + x */
      
      /* Copy fp into f */ 
      spartan_dcopy(fp,f,n);  /* f <- fp */
      obj = obj_true;         /* obj_true is the new obj */ 

      /* Compute the Jacobian */ 
      if (memmgmt) J = spartan_spfree(J);
      status = usrjac(&J,x,n);
      if(status){
        /* Clean up */ 
        spartan_free(work);
        return SPARTAN_USER_ERROR;
      }
      
      /* Compute a scaling matrix D */
      if(scaling) spartan_norms(J,d,dmin,dmax);
      
      /* Compute the gradient and inf norm  */ 
      spartan_zero(grad,n);       /* grad <- 0  */
      spartan_gatxpy(J,f,grad);   /* grad = J'f */
      finf = spartan_norminf(f,n); 

      /* Test for convergence to a local minmum */ 
      if(spartan_dnrm2(grad,n) <= jacstol && finf >= sqrt(tol)){
        converged = 1;
        status    = SPARTAN_LOCAL_MIN;
        snprintf(output,sizeof(char)*1024,"Warning: Converged to local minimum of merit function");
        spartan_print(output);
      }
      if(finf <= tol){
        converged = 1; 
        status    = SPARTAN_OK;
      }
      
    }
    k++;

    /* Print status information */
    snprintf(output,sizeof(char)*1024,"%3d %8.1e %8.1e %8.1e %8.1e %8.1e", \
             (int) k,finf,rho,normDp,Delta_old,fabs(obj_true - obj_model));
    spartan_print(output);
  }
  if (k==itermax && !converged) status = SPARTAN_ITER_LIMIT;
  
  /* Clean up */
  spartan_free(work);
  if (memmgmt) spartan_spfree(J);
  return(status);
}


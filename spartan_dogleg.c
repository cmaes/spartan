#include "spartan.h"
#include <stdio.h>

/* p = spartan_dogleg(f,J,Delta,d)

Computes an approximate solution to the trust-region subproblem 
 minimize    1/2 || J*p + f ||_2^2 
   p 
 subject to   || D*p|| <= Delta 

where D = diag(d) 
*/

#define CAUCHY_ON_BOUNDARY 1
#define NEWTON_IN_REGION   2
#define JACOBIAN_SINGULAR  3
#define DOGLEG_STEP_TAKEN  4
int spartan_dogleg(const double *f, const spartan_sparse *J, const double Delta,
                   const double *d, double *p){

  spartan_index n;
  double *g, *w;
  double eta, rho, lambda, lambdau, lambdac;

  n = J->n;

  /*  Note this problem is equivalent to 
      minimize    psi(p) = p'*g + 1/2*p'*H*p 
          p 
      subject to  || D*p || <= Delta 
      
      where the gradient g = J'*f and the 
      Hessian H = J'*J 
  */ 

  /* Allocate memory for g and w */
  g = spartan_calloc(n,sizeof(double));
  w = spartan_calloc(n,sizeof(double));

  /* We first compute the Cauchy point pc. This is the point that solves
     the problem
     minimize      psi(pc)
        pc 
     subject to    || D*pc || <= Delta 
                   pc  = -lambda*g 
      For some scalar lambda >= 0 
  */

  /* Compute the gradient */ 
  spartan_gatxpy(J,f,g);  /* g = J'*f */ 

  /* Compute rho = g'*H*g */
  spartan_gaxpy(J,g,w);   /* w = J*g */
  rho = spartan_ddot(w,w,n); /* rho = w'*w */

  /* Compute eta = || D*g || */
  spartan_dmult(d,g,w,n);    /* w = D*g */
  eta = spartan_dnrm2(w,n); 

  lambdau = spartan_ddot(g,g,n)/rho; /* lambdau = g'*g/rho */ 
  
  lambdac = Delta/eta;

  if (lambdau*eta >= Delta){
    spartan_zero(p,n);               /* p <- 0 */
    spartan_daxpy(-lambdac,g,p,n);   /* p = -lambdac*g */
    
    /* Clean up */ 
    spartan_free(w);
    spartan_free(g);

    /* spartan_print("Cauchy on boundary"); */
    return(CAUCHY_ON_BOUNDARY);
  }
  else {
    double *pc,*b;
    int solve_status;

    pc = spartan_calloc(n,sizeof(double));
    spartan_daxpy(-lambdau,g,pc,n);  /* pc = -lambdau*g */

    b  = spartan_calloc(n,sizeof(double));

    /* Compute the Newton step */
    spartan_daxpy(-1.0,f,b,n);               /* b = -f */
    solve_status = spartan_linsolve(J,b,p);  /*J*p = -f */

    if (solve_status == SPARTAN_SINGULAR_JACOBIAN){
      /* The Jacobian was singular take the Cauchy step */
      spartan_dcopy(pc,p,n);   /* p <- pc */ 

      /* Clean up */ 
      spartan_free(b);
      spartan_free(pc);
      spartan_free(w);
      spartan_free(g);

      /* spartan_print("Jacobian singular"); */
      return (JACOBIAN_SINGULAR);
    }
    else{
      spartan_dmult(d,p,w,n);   /* w = D*p */
      if (spartan_dnrm2(w,n) <= Delta){
        /* The Newton step is inside the trust-region. Take it.*/ 
        
        /* Relax. p already contains the Newton step */ 

        /* Clean up */
        spartan_free(b);
        spartan_free(pc);
        spartan_free(w);
        spartan_free(g);

        /* spartan_print("Netwon step taken"); */
        return(NEWTON_IN_REGION);
      }
      else{
        /* Compute the dogleg step */
        double alpha, beta, gamma, tau; 
        double *u; 

        /* Set w = pn - pc */ 
        spartan_dcopy(p,w,n);         /* w <- pn */ 
        spartan_daxpy(-1.0,pc,w,n);   /* w <- -1.0*pc + pn */ 
        
        /* Compute alpha = || D*(pn - pc)||^2 */ 
        spartan_dmult(d,w,b,n);         /* b <- d.*(pn - pc) */ 
        alpha = spartan_ddot(b,b,n);    /* alpha = b'*b */ 

        /* Compute beta = (D*pc)'*D*(pn - pc) */ 
        u = g;  /* Since g is no longer needed reuse its memory */ 
        spartan_dmult(d,pc,u,n);        /* u <- d.*pc */ 
        beta = spartan_ddot(u,b,n);     /* beta = u'*b */
        
        /* gamma = Delta^2 - || D*pc ||^2 */ 
        gamma = Delta*Delta - lambdau*lambdau*eta*eta;

        if (beta <= 0)
          tau = -beta + sqrt(beta*beta + alpha*gamma)/alpha;
        else
          tau = gamma/(beta + sqrt(beta*beta + alpha*gamma));
       
        /* Set p to the dogleg step */ 
        spartan_daxpy(tau,w,pc,n);     /* pc <- tau*(pn - pc) + pc */
        spartan_dcopy(pc,p,n);         /* p  <- pc + tau*(pn - pc) */
        

        /* Clean up */ 
        spartan_free(b);
        spartan_free(pc);
        spartan_free(w);
        spartan_free(g);
        
        {
          char buf[1024];
          sprintf(buf,"Dogleg step taken tau:%e alpha:%e beta:%e gamma:%e",tau,alpha,beta,gamma);
          /* spartan_print(buf); */
        }
        return(DOGLEG_STEP_TAKEN);
      }
    }

  }
  
}

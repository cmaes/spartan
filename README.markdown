SPARTAN: A Sparse Trust-Region Algorithm for Nonlinear Equations 
----------------------------------------------------------------

SPARTAN solves square systems of nonlinear equations, F(x) = 0, with
sparse derivatives. It is written in ANSI/ISO C, with a Scilab and
Matlab interface, and released under the GPL. SPARTAN relies on Tim
Davis’s UMFPACK routines to efficiently solve a sequence of sparse
linear systems, and Coleman, Garbow, and More's DSM subroutine to
compute a column coloring of the Jacobian matrix.

More details about Spartan are available in the 
[documentation](https://github.com/cmaes/spartan/blob/master/doc/spartan.pdf?raw=true)

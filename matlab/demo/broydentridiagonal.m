function [f,J] = broydentridiagonal(x)
% Broyden's tridiagonal function is:
% f_1(x) = (3 - h x_1)x_1 - 2 x_2 + 1, 
% f_i(x) = (3 - h x_i)x_i - x_{i-1} - 2 x_{i+1} +1 
%      i = 2, ..., n-1 
% f_n(x) = (3 - h x_n)x_n - x_{n-1} + 1
% h is a parameter h = 0.5 or h = 2
% x0 = (-1, ..., -1)^T 

n = size(x,1);
h = 0.5; 

f = zeros(n,1);
f(1) = (3 - h*x(1))*x(1) - 2*x(2) + 1;
for i=2:n-1
  f(i) = (3 - h*x(i))*x(i) - x(i-1) - 2*x(i+1) + 1;
end
f(n) = (3 - h*x(n))*x(n) - x(n-1) + 1;

if nargout > 1 
  nnzj = 2 + 3*(n-2) + 2; 
  irow = zeros(nnzj,1);
  jcol = zeros(nnzj,1);
  vals = zeros(nnzj,1);
  
  k = 1; 
  for i=1:n 
    if (i > 1)
      irow(k) = i;
      jcol(k) = i-1;
      vals(k) = -1;
      k = k+1; 
    end
    
    irow(k) = i;
    jcol(k) = i; 
    vals(k) = 3 - 2*h*x(i);
    k = k+1; 

    if (i < n)
      irow(k) = i; 
      jcol(k) = i+1;
      vals(k) = -2;
      k = k+1;
    end
  end
  
  % form a matlab sparse matrix 
   J = sparse(irow,jcol,vals,n,n);

  % form a scilab sparse matrix 
  %J = sparse([irow jcol],vals,[n,n]);
end 
 
 

function [F,J] = generalrosenbrock(x)
% Derived from the generalized n-dimensional 
% Rosenbrock function, returns a dense Jacobian 
n = length(x);         
if n == 0, error('Input vector, x, is empty.'); end
if mod(n,2) ~= 0, 
   error('Input vector, x ,must have an even number of components.');
end
%Evaluate the vector function
odds  = 1:2:n;
evens = 2:2:n;
F = zeros(n,1);
F(odds,1)  = 1-x(odds);
F(evens,1) = 10.*(x(evens)-x(odds).^2); 
%Evaluate the Jacobian matrix if nargout > 1
if nargout > 1
   c = -ones(n/2,1);    C = sparse(odds, odds,c,n,n);
   d = 10*ones(n/2,1);  D = sparse(evens, evens,d,n,n);
   e = -20.*x(odds);    E = sparse(evens, odds,e,n,n);
   J = C + D + E;
end

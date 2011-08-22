function [f,J] = powellsingular(x)
    n = size(x,1);
    if n ~= 2, error('bad size of x'); end

    f = zeros(n,1);
    f(1) = x(1);
    f(2) = 10*x(1)/(x(1) + 0.1) + 2*x(2)^2; 

    if nargout > 1
        J = [ 1                      0; 
              1/((x(1)+0.1)^2)  4*x(2)];
        J = sparse(J);
    end

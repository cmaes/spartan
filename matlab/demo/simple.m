function [f,J] = simple(x) 
       f = (x-1).^2;
       if nargout > 1 
           J = diag(sparse(2.*(x-1)));
       end
function [f,J] = singularbroyden(x)
   n = size(x,1);
   h = 0.5;

   f = zeros(n,1);

   % this is the broyden function 
   f(1) = ((3-h*x(1))*x(1) - 2*x(2) + 1);
   for i=2:n-1
     f(i) = ((3-h*x(i))*x(i) - x(i-1) - 2*x(i+1) + 1);
   end
   f(n) = ((3 - h*x(n))*x(n) - x(n-1) + 1);
   
   % however, we need to square each f_i 
   % we will do this after we compute the jacobian


   % if nargout > 1 compute the sparse Jacobian
   if nargout > 1 
     % the jacobian is tridiagonal 
     nnzjac = 2 + 2 + 3*(n-2);
     
     irow = zeros(nnzjac,1);
     jcol = zeros(nnzjac,1);
     vals = zeros(nnzjac,1);
     k = 1; 
     
     for i=1:n
       % \partial f_i / \partial x_{i-1}
       if i > 1 
	 irow(k) = i;
	 jcol(k) = i-1;
	 vals(k) = 2*f(i)*(-1);
	 k = k + 1; 
       end
       
       % \partial f_i / \partial x_i 
       irow(k) = i;
       jcol(k) = i;
       vals(k) = 2*f(i)*(3 -2*h*x(i));
       k = k + 1;

       % \partial f_i / \partial x_{i+1}
       if i < n
	 irow(k) = i;
	 jcol(k) = i+1;
	 vals(k) = 2*f(i)*(-2);
	 k = k + 1;
       end
     end
     % create a scilab sparse matrix 
     
     J = sparse(irow,jcol,vals,n,n);
   end

   % square f now
   f = f.^2;
end
 

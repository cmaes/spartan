function [p,status] = spartan_dogleg(f,J,Delta,d) 
% p = spartan_dogleg(f,J,Delta,d) 
% Computes an approximate solution to the trust-region subproblem 
% minimize       1/2 || f + J*p ||_2^2  
%    p 
% subject to      || D*p || <= Delta 
% where D = diag(d)
% using the dogleg method

% Algorithm Parameters 
singulartol = 1e-6;

% Check some sizes 
[m,n] = size(J);
if m ~= n, error('Jacobian must be square'); end

% Note this problem is equivalent to 
% minimize    psi(p) = p'*g + 1/2*p'*H*p 
%    p 
% subject to  || D*p || <= Delta 
%
% where the gradient g = J'*f and the 
% Hessian H = J'*J
% 

% We first compute the Cauchy point pc. This is the point that solves
% the problem
% minimize      psi(pc)
%    pc 
% subject to    || D*pc || <= Delta 
%               pc  = -lambda*g 
% For some scalar lambda >= 0 

% Compute the gradient 
g = J'*f;

% rho = g'*H*g
w = J*g;  rho = w'*w; 

eta = norm(d.*g); % eta = || D*g||

% lambdau is the step to the unconstrained minimizer
% of psi(pc) = -lambda g'*g + lambda^2/2 g'*H*g
lambdau = g'*g/rho;

% lambdac is the step to the trust-region boundary 
lambdac = Delta/eta;

if lambdau*eta >= Delta
    % The unconstrained minimizer lies outside of the trust-region. Choose
    % lambda so that the Cauchy point lies on the boundary of the
    % trust-region. In this case the Cauchy point is a good
    % approximate solution.
    p = -lambdac*g;    
    status = 'Cauchy hit boundary';
else 
    % The unconstrained minimizer lies inside the trust-region. 
    pc = -lambdau*g;

    % Try to calculate the Newton step pn that satisfies 
    % H*pn = -g or (J'*J)*pn = -(J'*f) 
    % Note that the solution pn solves the least-squares problem 
    % minimize || J*pn + f ||
    %   pn 

    % Disable singular matrix warnings 
    s = warning('off','MATLAB:singularMatrix');
    
    % Compute the solution to the least-squares problem, via
    % Matlab's sparse QR 
    %perm = colamd(J);              % Compute a fill reducing column ordering 
    %[w,R] = qr(J(:,perm),-f);      % J*Dinv*P = Q*R and w = Q'*(-f)  
    %pn = R\w;                      % pn is soln to min || J*Dinv*P*x + f ||
    %pn = pn(perm);                 % Compute pn = P*x

    % Compute the minimum 2-norm soln to the rank-deficent least-squares 
    % problem using Tim Davis SuiteSparse QR 
    % [c,R,perm] = spqr(J,-f,struct('econ',0,'permutation','vector'));
    %if size(R,1) ~= n 
    %  fprintf('min2norm\n')
    %  [y] = spqr_solve(R,c,struct('solution','min2norm'));
    %  pn = y(perm);
    %else 
    % fprintf('exact %d %d\n',size(perm,1),size(perm,2));
    % y = R\c;
    % y(perm) = y;
    % pn = y;
    %end
   
    % Compute the Newton step via Matlab's sparse LU. 
    % If J is singular pn will have inf or nan
    pn = J\(-f);
   
    % Restore warning state
    warning(s);
    
    % If J is nonsigular then the least-squares problem will have zero
    % residual
    residual = norm(J*pn + f);
    %fprintf('residual:%e ',residual);

    % If the residual is large then J was singular. In this case we take
    % the solution to be the Cauchy point.
    if   any(isinf(pn)) || any(isnan(pn)) % || residual > singulartol
        p = pc;                % This is the scaled solution 
        status = sprintf('Jacobian singular: %e',residual);
    else
        if norm(d.*pn) <= Delta 
            % The Newton step is inside the trust-region radius.  So it solves the
            % trust-region problem.
            p = pn;            % This is the scaled solution
            status = sprintf('Newton step taken %e %e',residual,norm(pn));
        else 
           % The Netwon step is outside the trust-region radius.  We compute the
           % dogleg step pd that satisifes
           % maximize      tau 
           % subject to    pd = pu + tau*(pn - pu)
           %               || D*pd || <= Delta
           %               0 <= tau <= 1 
           % The equation || D*pd || = Delta^2 yields the 
           % a quadratic equation in tau, which we solve .

           % pu is the Cauchy point since ||D*pc|| <= Delta 
           pu = pc; 

           % These steps solve the quadratic equation. 
           % omega(tau) = || D*pd(tau) ||^2 - Delta^2 = 0 
           % omega(tau) = alpha*tau^2 + 2*beta*tau + ||D*pu||^2 - Delta^2 = 0 
           % where alpha = || D*(pn - pu) ||^2 >= 0 and
           %       beta  = (D*pu)'*(D*(pn - pu))
           % let   gamma = Delta^2 - || D*pu||^2 
           % 
           % Note that omega(0) = || D*pu ||^2 - Delta^2 < 0 
           % and  that omega(1) = || D*pn ||^2 - Delta^2 > 0 
           % so a root of omega lies in the interval (0,1)
           % 
           % Also, note that -alpha*(||D*pu||^2 - Delta^2) > 0
           % so both roots are real 

           w   = pn - pu;
           Dw  = d.*w; 
           Dpu = d.*pu;
           
           alpha = Dw'*Dw;  % (pn - pu)'*D^2*(pn - pu)
           beta  = Dpu'*Dw; % pu'*D^2*(pn-pu) 
           gamma = Delta^2 - (lambdau*eta)^2; % Delta^2 - || D*pu||^2

           if beta <= 0 
               % This is the (-b + sqrt(b^2 - 4ac))/2a case. 
               % Since a = alpha >= 0 and b < 0 the root is
               % positive. 
               
               % We compute the root this way to avoid subtraction,
               % the factor of two has also been removed 
               tau = -beta + sqrt(beta^2 + alpha*gamma)/alpha;
           else 
               % This is the (-b - sqrt(b^2 - 4ac))/2a case. 
               tau = gamma/(beta + sqrt(beta^2 + alpha*gamma));
           end

           % Set p to the dogleg step 
           p = pu + tau*(pn - pu);   
           status = sprintf('Dogleg taken res:%e norm:%e Delta:%e tau:%e diff:%e alpha:%e beta:%e gamma:%e',...
                            residual,norm(d.*p),Delta,tau,Delta-norm(d.*p),alpha,beta,gamma);
        end
    end
end 
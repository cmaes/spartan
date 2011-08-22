function [x] = spartan(func,x) 
% [x] = spartan(func,x0)
% Solves the nonlinear system f(x) = 0 using a trust-region algorithm

% We seek to minimize the merit function 
% phi(x) = 1/2 || f(x) ||_2^2 

n = size(x,1);

% Parameters for the algorithm 
DeltaMax = 100;
Delta    = 1e+0;
eta      = 1e-3; 
tol      = 1e-6; 
jacstol  = 1e-12;
itmax    = 100; 
dmin     = eps;
dmax     = 1/eps;

k = 0; 
converged = 0;
fprintf('itn   ||f||    rho     ||Dp||    Delta    err\n');

% Evaluate the function and Jacobian at the starting point 
[f,J] = func(x);
d = ones(n,1);
% Compute scaling for the Jacobian 
for j=1:n 
    d(j) = norm(J(:,j));
end
d = max(d,dmin);
d = min(d,dmax);
d = ones(n,1);

fprintf('%3d %8.1e                   %8.1e\n',k,norm(f,inf),Delta);
while ~converged && k < itmax
    
    % Calculate p as an approximate solution to the local quadratic
    % approximation to phi, subject to a trust-region constraint
    %
    % minimize psi(p) = 1/2 || f + J p ||_2^2 
    %    p 
    % subject to    || D*p ||_2 <= Delta 
    %
    % where f = f(x_k) and J is the Jacobian f'(x_k)
    %
    % We solve this problem using the dogleg method 

    
    p = spartan_dogleg(f,J,Delta,d);
    normDp  = norm(d.*p);
    Deltao = Delta;
    str = '';

    % Compute the function value predicted at this step 
    fpred = f + J*p;

    % Compute the function value at this step
    fp = func(x + p);

    % Evaluate rho the ratio of the actual to predicted reduction
    rho = (f'*f - fp'*fp) / ( f'*f - fpred'*fpred );
    
    if rho < 0.25
        % The models prediction was not good and we did not acheive a good
        % reduction. Decrease the trust-region radius
        Delta = 0.25*Delta; 
    else

        if rho > 0.75 && abs(normDp - Delta) < 1e-8 
            % There is good agreement between the model and the
            % function (rho is close to 1), and the step is hitting
            % the trust-region bound. Thus, it is safe to increase
            % the trust-region radius.
            Delta = min(2*Delta, DeltaMax);
        end

    end 

    if rho > eta 
        % We acheived a significant reduction. 
        % Take the step 
        x = x + p; 
        
        % Update the function and the Jacobian
        [f, J] = func(x);

        % Compute scaling for the Jacobian
        for j=1:n 
            d(j) = norm(J(:,j));
        end
        d = max(d,dmin);
        d = min(d,dmax);
        d = ones(n,1);
        
        % Test for convergence 
        if norm(f,inf) <= tol 
            converged = 1;
        end
        
        % Test for convergence to local minimum 
        if norm(J'*f) <= jacstol && norm(f,inf) >= 100*tol
            converged = 1;
            warning(['Converged to a local minimum of merit function ' ...
                     'but not a solution to system of equation']);
        end
    end
    % Otherwise we compute a new p with a new radius. 

    k = k + 1; 
    % Print status information 
    fprintf('%3d %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %8.1e %s\n',...
            k,norm(f),rho,normDp,Deltao,norm(fpred - fp),min(d),max(d),str);
end 

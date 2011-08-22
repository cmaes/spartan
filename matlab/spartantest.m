n = 2; 
x0 = zeros(n,1);
fprintf('\nSolving simple diagonal Jacobian\n');
f = @(x) simple(x);

tic; 
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x-1) = %e\n',norm(x-1));
fprintf('time      = %g\n',timetaken);

n = 1000; 
x0 = -ones(n,1);
fprintf('\nSolving function6 %d\n',n);
f = @(x) function6(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x)   = %e\n',norm(x));
fprintf('time      = %g\n',timetaken);

n = 1000; 
x0 = -ones(n,1);
fprintf('\nSolving trigexp %d\n',n);
f = @(x) trigexp(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x)   = %e\n',norm(x));
fprintf('time      = %g\n',timetaken);

n = 1000;
x0 = -ones(n,1);
fprintf('\nSolving Broyden tridiagonal %d\n',n);
f = @(x) broydentridiagonal(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x)   = %e\n',norm(x));
fprintf('time      = %g\n',timetaken);


n = 1000;
x0 = -ones(n,1);
fprintf('\nSolving Singular Broyden tridiagonal %d\n',n);
f = @(x) singularbroyden(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x)   = %e\n',norm(x));
fprintf('time      = %g\n',timetaken);

n = 2^10;
x0 = zeros(n,1);
x0(1:n,1) = -1.9;
x0(2:2:n,1) = 2; 

fprintf('\nSolving Genralized Rosenbrock %d\n',n);
f = @(x) generalrosenbrock(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x-1) = %e\n',norm(x-1));
fprintf('time      = %g\n',timetaken);

n = 2;
x0 = [1.8016; 0];

fprintf('\nSolving Powell''s singular function %d\n',n);
f = @(x) powellsingular(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x)   = %e\n',norm(x));
fprintf('time      = %g\n',timetaken);


n = 1000;
x0 = -ones(n,1);

fprintf('\nSolving nlsf1 %d\n',n);
f = @(x) nlsf1(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x)   = %e\n',norm(x));
fprintf('time      = %g\n',timetaken);

n = 1;
x0 = 0*ones(n,1);

fprintf('\nSolving minpackatol %d\n',n);
f = @(x) minpackatol(x);

tic
x = spartan(f,x0);
timetaken = toc;
fprintf('norm(f)   = %7.1e\n',norm(f(x)));
fprintf('norm(x)   = %e\n',norm(x));
fprintf('time      = %g\n',timetaken);
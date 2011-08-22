function [f,J] = minpackatol(x)
    if all(size(x) ~= [1,1]) error('bad x dimension'), end
    f = (3*x - 10)*exp(10*x);
    J = 3*exp(10*x) + (3*x - 10)*exp(10*x)*10;
    J = sparse(J);
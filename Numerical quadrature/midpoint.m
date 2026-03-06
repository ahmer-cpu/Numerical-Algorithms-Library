function [int,error] = midpoint(func,mesh,H,n,true_int)
%MIDPOINT Summary of this function goes here
%   Detailed explanation goes here
h = H/2.0;
c = mesh(1:n) + h;          % Broadcasting for c_i values
fvals = func(c);

int = H*sum(fvals);

error = true_int-int;
end


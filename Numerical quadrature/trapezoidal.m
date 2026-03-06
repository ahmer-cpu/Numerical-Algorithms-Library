function [int,error] = trapezoidal(func,mesh,H,n,true_int)
%TRAPEZOIDAL Summary of this function goes here
%   Detailed explanation goes here
fvals = func(mesh);

int = 2.0*sum(fvals(2:n));
int = int+fvals(1)+fvals(n+1);
int = int*(H/2.0);

error = true_int-int;
end


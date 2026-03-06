function [int,error] = opentwopoint(func,mesh,H,n,true_int)
%OPENTWOPOINT Summary of this function goes here
%   Detailed explanation goes here
h = H/3.0;
xprime = mesh(1:n) + h;         % Broadcasting for x_i' values
xprimeprime = mesh(2:n+1) - h;  % Broadcasting for x_i'' values
fvals_left = func(xprime);
fvals_right = func(xprimeprime);

int = sum(fvals_left)+sum(fvals_right);
int = int*(H/2.0);

error = true_int-int;
end


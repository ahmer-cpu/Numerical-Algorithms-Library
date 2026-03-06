function [int,error] = gausslege(func,mesh,H,n,true_int)
%OPENTWOPOINT Summary of this function goes here
%   Detailed explanation goes here

H_left = (H/2.0)*(1-(1/sqrt(3)));
H_right = (H/2.0)*(1+(1/sqrt(3)));

xprime = mesh(1:n) + H_left;            % Broadcasting for x_i' values
xprimeprime = mesh(1:n) + H_right;      % Broadcasting for x_i'' values

fvals_left = func(xprime);
fvals_right = func(xprimeprime);

int = (H/2) * sum(fvals_left + fvals_right);

error = true_int-int;
end
function [int,error] = simpson(func,mesh,H,n,true_int)
%SIMPSON Summary of this function goes here
%   Detailed explanation goes here
fvals = func(mesh);
h = H/2.0;
c = mesh(1:n) + h;          % Broadcasting for c_i values
fc_vals = func(c);

int = 2.0*sum(fvals(2:n));
int = int+(4.0*sum(fc_vals));
int = int+fvals(1)+fvals(n+1);
int = int*(H/6.0);

error = true_int-int;
end


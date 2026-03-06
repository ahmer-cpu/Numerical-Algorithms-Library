function [int,error] = leftrectangle(func,mesh,H,n,true_int)
%LEFTRECTANGLE Summary of this function goes here
%   Detailed explanation goes here
fvals = func(mesh(1:n));

int = sum(H*fvals);

error = true_int-int;
end


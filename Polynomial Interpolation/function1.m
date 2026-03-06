function fvals = function1(xvals)
%FUNCTION1 Outputs the values of function 1 using the product method
%   This routine takes in the parameters
%   p and d to evaluate (x-p)**d for all xvals
d = 9;
p = 2;
nodes = zeros(1, d);
nodes(:) = p;         
fvals = polyProduct(1,nodes,xvals);
end


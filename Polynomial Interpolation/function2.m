function fvals = function2(xvals)
%FUNCTION2 Outputs the values of function 2 using the product method
%   This routine takes in the parameters
%   d to evaluate \prod(1,d) (x-i) for all xvals

d = 10;
nodes = 1:d;
fvals = polyProduct(1,nodes,xvals);
end


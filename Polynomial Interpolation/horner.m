function p = horner(alphas,xnodes,xvals)
%HORNER Evaluate the polynomial in the Newton Basis form
%   This functions takes in as input the alphas and xnodes that define a
%   polynomial in the newton basis and evaluates it at the points in xvals
%   and returns the matrix/vector of evaluations

p = zeros(size(xvals));
p(:) = alphas(end);
for i = length(alphas)-1 : -1: 1
    p = (p .* (xvals - xnodes(i))) + alphas(i);
end
end
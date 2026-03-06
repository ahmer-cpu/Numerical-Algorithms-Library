function p = polyProduct(alpha,nodes,xvals)
%POLYPRODUCT Polynomial Evaluation for the product form
%   This function takes in paramaters defining a polynomial in product
%   form, alpha and the nodes and evaluates it on the array xvals, returning an array of evaluations

p = zeros(size(xvals));
p(:) = alpha;
for i = 1:length(nodes)
    p = p .* (xvals - nodes(i));
end



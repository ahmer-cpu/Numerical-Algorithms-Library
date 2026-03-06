function fvals = function3(xvals,nodes)
%FUNCTION3 Outputs the values of function 3
%   This function is the last lagrange basis for a set of nodes (parameter.
nodes = nodes(1:end-1);

fvals = polyProduct(1,nodes,xvals);

div = polyProduct(1,nodes,nodes(end));

fvals = fvals / div;
end


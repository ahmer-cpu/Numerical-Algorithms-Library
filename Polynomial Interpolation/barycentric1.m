function [gamma,fvals] = barycentric1(func, n, nodes)
%BARYCENTRIC1 Computes the array of Barycentric coefficients gammas and
%the function values at the nodes
%   This routine takes in as parameters a function, the number of nodes n,
%   and the array of nodes and computes the coefficients gammas required
%   for interpolation in Barycentric form 1, as well as the function values
%   at the nodes

n = single(n);
nodes = single(nodes);
fvals = single(func(single(nodes)));

if length(nodes) < n
    error('nodes must have at least n elements.');
end

if n == 1
    m = 1;  % Only one node: trivial case.
    return;
end

m = single(zeros(1, n));

m(1) = nodes(1) - nodes(2);
m(2) = nodes(2) - nodes(1);

for k = 2:(n-1)
    p = single(1);
    for i = 0:(k-1)
        t = nodes(i+1) - nodes(k+1);
        m(i+1) = t * m(i+1);
        p = -t * p;
    end
    m(k+1) = p;
end

gamma = single(1) ./ single(m);
end
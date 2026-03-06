function [dd,fvals] = newtondd(func,nodes)
%NEWTONDD Compute divided difference table
%   The function takes as input a function parameter and x values and
%   outputs the top row/left column of the divided differences table for
%   interpolation assuming the x values are indexed from top to bottom in the array

n = length(nodes);
n = single(n);

fvals = func(nodes);
fvals = single(fvals);

nodes = single(nodes);

dd = fvals;
for i = 1:n-1
    dd(i+1:end) = dd(i+1:end) - dd(i:end-1);
    D = nodes(i+1:end) - nodes(1:end-i);
    dd(i+1:end) = (dd(i+1:end)) ./ D;
end

end


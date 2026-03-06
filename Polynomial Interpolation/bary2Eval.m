function p = bary2Eval(xvals,fvals,nodes,betas,n)
%BARY2EVAL Evaluates polynomial in Barycentric 2 form
%   This routine takes in as parameters fvals, nodes, and betas and
%   evaluates the polynomial at an array of xvals

xvals = single(xvals);
fvals = single(fvals);
nodes = single(nodes);
betas = single(betas);
n = single(n);


numer = single(zeros(size(xvals)));
denom = single(zeros(size(xvals)));
for j = 1:n
    xdiff = xvals - nodes(j);
    temp = betas(j) ./ xdiff;
    numer = numer + (temp*fvals(j));
    denom = denom + temp;
end
p = numer ./ denom;

tol = 1e-12;

dists = abs( repmat(xvals, length(nodes), 1)' - repmat(nodes, length(xvals), 1) );

% For each xval (each row), find the minimum distance and its corresponding node index.
[minDist, minDistIdx] = min(dists, [], 2);

% Identify which xvals are too close to a node.
closeMask = (minDist < tol);

% For each xval that is too close, assign p at that xval the value from fvals 
% corresponding to the index of the closest node.
p(closeMask) = fvals(minDistIdx(closeMask));
end


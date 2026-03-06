function [beta,fvals] = barycentric2(func,n,flag)
%BARYCENTRIC2 Gives the betas and the function values computed at certain
%nodes for Barycentric 2 form.
%   This routine takes in as parameters a function, the number of nodes,
%   and a flag indicating the type of nodes and returns the betas and fvals
if nargin < 3
    flag = 'uniform';
end

switch lower(flag)
    case 'uniform'
        beta = single(zeros(1,n));
        beta(1) = single(1);
        for i = 0: (n-2)
            beta(i+2) = -beta(i+1) * ((n-i )/(i+1));
        end
        nodes = uniform(-1,1,n);
        
    case 'cheby1'
        j = 0:n-1;
        beta = ((-1).^j) .* sin((2*j + 1)*pi / (2*n + 2));
        nodes = chebyshev1(-1,1,n);
        
    case 'cheby2'
        beta = [0.5; ones(n-2,1); 0.5].*(-1).^((0:n-1)');
        nodes = chebyshev2(-1,1,n);
        
    otherwise
        error('Unknown ordering flag: %s. Use "uniform", "cheby1", or "cheby2".', flag);
end

fvals = single(func(nodes));

end

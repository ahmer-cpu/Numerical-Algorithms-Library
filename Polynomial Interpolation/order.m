function x_ordered = order(x,flag)
%ORDER Compute ordered array for the given ordering type
%   This function takes an as input an unordered array and a flag
%   indicating the type of ordering, and returns the ordered array (Leja,
%   descending, or ascending)

if nargin < 2
    flag = 'ascend';
end

switch lower(flag)
    case 'ascend'
        x_ordered = sort(x, 'ascend');
    case 'descend'
        x_ordered = sort(x, 'descend');
    case 'leja'
        x_ordered = lejaOrder(x);
    otherwise
        error('Unknown ordering flag: %s. Use "ascend", "descend", or "leja".', flag);
end
end

function leja = lejaOrder(x)
n = length(x);
leja = zeros(size(x));

[~, idx] = max(abs(x));
leja(1) = x(idx);
x(idx) = [];

for k = 2:n
    products = zeros(size(x));
    for j = 1:length(x)
        products(j) = prod(abs(x(j) - leja(1:k-1)));
    end
    [~, maxIdx] = max(products);
    leja(k) = x(maxIdx);
    x(maxIdx) = [];
end
end

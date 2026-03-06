function x = chebyshev1(a,b,n)
j = 0:n-1;
x = cos( (2*j + 1)*pi / (2*n + 2) );
x = lineartransform(x,-1,1,a,b);
end

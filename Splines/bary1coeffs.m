function [gammainvs] = bary1coeffs(xvals,n)

%-------------------------------------------------------
%
% This code computes the inverses of the coefficients
% for the barycentric form 1 (modified Lagrange) of
% the interpolating polynomial on the mesh given
% in the n+1 by 1 column vector xvals.
% They are returned in the n+1 by 1 column vector gammainvs.
% n is the degree of the interpolating polynomial to be formed. 
%
%-------------------------------------------------------


  np1=n+1;
gammainvs = zeros(np1,1); 


%-------------------------------------------------------
% initialize to the linear coefficients 
%-------------------------------------------------------

gammainvs(1) = (xvals(1) - xvals(2));  %   (x_0 - x_1)
gammainvs(2) = (xvals(2) - xvals(1));  %   (x_1 - x_0)
d=1;  % degree for whom the parameters are in gammainvs

%-------------------------------------------------------
% from degree 2 to n
% add effect of x_d, i.e., xvals(d+1)
% to gammainvs
%
%-------------------------------------------------------
for d=2:n
  p=1.0;
  for i= 0:d-1;
      t=(xvals(i+1) - xvals(d+1));
      gammainvs(i+1)=gammainvs(i+1)*t;  % m_i^d = t x m_i^{n-1}
      p=-p*t;
  end
  gammainvs(d+1)= p;
end

return
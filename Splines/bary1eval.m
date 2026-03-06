function [pvals,kappa_numer,lebesgue_vals] = bary1eval(xvals,yvals,gammainvs,n,xcheck,ncheck,kappa_numer,lebesgue_vals)

%-------------------------------------------------------
%
% This code evaluates the interpolating polynomial p(x)
% of degree d=n in barycentric form 1 (modified Lagrange) 
% for the data pairs in (xvals(i),yvals(i)) specifed by
% the barycentric form 1 coefficients in gammainvs
% computed by bary1coeffs.m.
% xvals, yvals and gammainvs are all assumed n+1 by 1
% column vectors.
%
% The output is an ncheck by 1 column vector pvals
% containing pvals(i)=p(xcheck(i))
%
% The numerator of kappa(x,n,f) is returned
% in kappa_numer(1:ncheck)
%
% kappa(x,n,1) values are returned in lebesgue_vals(1:ncheck)
%
%-------------------------------------------------------


%-------------------------------------------------------
%
% EVALUATION ON A FINE GRID OF POINTS
%
% protection made if x = x_i for some i
% see Berrut and Trefethen for brief discussion
% exploits 0/0 = NAN but does not stop
%
%-------------------------------------------------------

omega = ones(ncheck,1);
exact = zeros(ncheck,1);
onevals = omega;
np1=n+1;

for i=1:np1
  xdiff=(xcheck - onevals*xvals(i));
  exact(xdiff==0) = i;
  omega = omega.*xdiff;
end

%
% The numerator of kappa(x,n,f) is computed and returned
% in kappa_numer(1:ncheck)
% and pvals the polynomial using the original data.



pvals = zeros(size(xcheck));
kappa_numer = zeros(size(xcheck));
lebesgue_vals= zeros(size(xcheck));

for i=1:np1
  xdiff=(xcheck - onevals*xvals(i));
  pvals = pvals + ((onevals*(yvals(i)/gammainvs(i))).*(omega./xdiff));
  kappa_numer = kappa_numer + (abs(onevals*(yvals(i)/gammainvs(i))).*abs(omega./xdiff));
  lebesgue_vals = lebesgue_vals  + (abs(onevals*(1.0/gammainvs(i))).*abs(omega./xdiff));
end

jfixes = find(exact);
pvals(jfixes)=yvals(exact(jfixes));
kappa_numer(jfixes)=abs(yvals(exact(jfixes)));
lebesgue_vals(jfixes)=abs(yvals(exact(jfixes)));

pvals = pvals';
kappa_numer = kappa_numer';
lebesgue_vals = lebesgue_vals';
return

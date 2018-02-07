function [cpoly,dcpoly] = characteristic_polyfit(x)
% characteristic_polyfit - Generate characteristic polynomial coefficients
% and associated derivatives for the given set of points
%
% Syntax:  [cpoly,dcpoly] = characteristic_polyfit(x)
%
% Inputs:
%    x - Set of points to fit
%
% Outputs:
%    cpoly - Matrix of with row i containing the polynomial coefficients
%            of the ith characteristic polynomial.
%    dcpoly - Matrix of with row i containing the polynomial coefficients
%             of the ith characteristic polynomial derivative.
%
% Example:
%    >> [cpoly,dcpoly] = characteristic_polyfit([0 0.5 1]);
%    >> cpoly
% 
%    cpoly =
%        2.0000   -3.0000    1.0000
%       -4.0000    4.0000   -0.0000
%        2.0000   -1.0000         0
%
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 08-Oct-2017

%------------- BEGIN CODE --------------

% Check inputs
if (size(x,1) == 1)
elseif (size(x,2) == 1)
    x = x.';
else
    error('Input x must be a vector');
end

% Number of points
npts = length(x);

% Find polynomial coefficients through points
cpoly = zeros(npts);
for n = 1:npts
    bvec = zeros(1,npts);
    bvec(n) = 1;
    cpoly(n,:) = polyfit(x, bvec, npts-1);
end

% Take derivatives of fit polynomials
dcpoly = zeros(npts);
for n = 1:npts
    dcpoly(n,2:npts) = cpoly(n,1:npts-1) .* [npts-1:-1:1];
end

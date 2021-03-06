function [M1,M2,M3] = se_legfilter_matrices(npts, coeffs)
% se_legfilter_matrices - Generate matrices for the Legendre mode filter.
%
% Syntax:  [M1,M2,M3] = se_legfilter_matrices(npts)
%
% Inputs:
%    npts - Order of the spectral element method
%    coeffs - Vector of length npts containing the filter coefficients
%
% Outputs:
%    M1 - Matrix of left coefficients (npts-1 x npts-1)
%    M2 - Matrix of center coefficients (npts-1 x npts-1)
%    M3 - Matrix of right coefficients (npts-1 x npts-1)
%
% Example: 
%    >> [L1,L2,L3] = se_legfilter_matrices(4);
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 19-Apr-2018

%------------- BEGIN CODE --------------

% Check arguments
if (~isscalar(npts))
    error('npts argument must be a scalar');
end
if (~isvector(coeffs))
    error('coeffs argument must be a vector');
end
if (length(coeffs) ~= npts)
    error('coeffs vector must have length npts');
end

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[myroots,myweights] = lglnodes(npts,0,1);

% Calculate modal transform matrix
D = zeros(npts,npts);
for n = 1:npts
    %D(:,n) = exp(1i * (n-1) * pi * myroots);
    D(:,n) = legendreP((n-1), (2*myroots-1));
end

% Build filter
MD = real(D * diag(coeffs) * inv(D));

MD(1,2:npts) = 0.5 * MD(1,2:npts);
MD(npts,1:npts-1) = 0.5 * MD(npts,1:npts-1);

% Exract submatrices
M1 = zeros(npts-1,npts-1);
M1(1,:) = MD(npts,1:npts-1);

M2 = MD(1:npts-1,1:npts-1);

M3 = zeros(npts-1,npts-1);
M3(:,1) = MD(1:npts-1,npts);

myweights;
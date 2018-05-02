function [M1,M2,M3] = se_bvfilter_matrices(npts, mu, porder, s)
% se_laplacian_matrices - Generate matrices for the Boyd-Vandeven filter.
%
% Syntax:  [M1,M2,M3] = se_bvfilter_matrices(npts)
%
% Inputs:
%    npts - Order of the spectral element method
%    mu - Non-dimensional filter coefficient
%    porder - Order of the Boyd-Vandeven filter
%    s - Filter lag
%
% Outputs:
%    M1 - Matrix of left coefficients (npts-1 x npts-1)
%    M2 - Matrix of center coefficients (npts-1 x npts-1)
%    M3 - Matrix of right coefficients (npts-1 x npts-1)
%
% Example: 
%    >> [L1,L2,L3] = se_bvfilter_matrices(4);
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

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[myroots,myweights] = lglnodes(npts,0,1);

% Calculate modal transform matrix
D = zeros(npts,npts);
for n = 1:npts
    D(:,n) = legendreP((n-1), (2*myroots-1));
end

% Coefficient weights at each node
coeffs = zeros(npts);
for k = 1:npts
    kx = k-1;
    if (kx < s)
        coeffs(k,k) = 1;
    else
        Omega = abs((kx-s)/(npts-s)) - 0.5;
        coeffs(k,k) = 0.5 * (1 - erf(2 * sqrt(porder) * sign(Omega) * sqrt(-log(1 - 4 * Omega^2) / 4)));
    end
end

% Build filter
MD = (1 - mu) * eye(npts) + mu * real(D * coeffs * inv(D));

MD(1,2:npts) = 0.5 * MD(1,2:npts);
MD(npts,1:npts-1) = 0.5 * MD(npts,1:npts-1);

% Exract submatrices
M1 = zeros(npts-1,npts-1);
M1(1,:) = MD(npts,1:npts-1);

M2 = MD(1:npts-1,1:npts-1);

M3 = zeros(npts-1,npts-1);
M3(:,1) = MD(1:npts-1,npts);

myweights;
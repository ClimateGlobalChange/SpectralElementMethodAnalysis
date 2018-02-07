function [M1,M2,M3] = se_laplacian_matrices(npts)
% se_laplacian_matrices - Generate matrices for the spectral element
% second derivative operator with unit element width.
%
% Syntax:  [M1,M2,M3] = se_laplacian_matrices(npts)
%
% Inputs:
%    npts - Order of the spectral element method
%
% Outputs:
%    M1 - Matrix of left coefficients (npts-1 x npts-1)
%    M2 - Matrix of center coefficients (npts-1 x npts-1)
%    M3 - Matrix of right coefficients (npts-1 x npts-1)
%
% Example (derivative of x^2): 
%    >> [L1,L2,L3] = se_laplacian_matrices(4);
%    >> [pts,~] = lglnodes(4, 0, 1);
%    >> L1*(pts(1:3)-1).^2+L2*pts(1:3).^2+L3*(pts(1:3)+1).^2
%
%    ans =
%       2.0000
%       2.0000
%       2.0000
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 08-Oct-2017

%------------- BEGIN CODE --------------

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[myroots,myweights] = lglnodes(npts,0,1);

% Find coefficients of characteristic polynomials
[~,difflegpoly] = characteristic_polyfit(myroots);

% Obtain mass matrix
M = myweights;
M(1) = 2 * M(1);
M(npts) = 2 * M(npts);

% Obtain differentiation matrix
D = zeros(npts,npts);
for n = 1:npts
for m = 1:npts
    for s = 1:npts
        D(n,m) = D(n,m) - myweights(s) ...
            * polyval(difflegpoly(n,:), myroots(s)) ...
            * polyval(difflegpoly(m,:), myroots(s));
    end
end
end

% Mirror integrals along edges
D(1,1) = 2 * D(1,1);
D(npts,npts) = 2 * D(npts*npts);

% Obtain product
MD = inv(diag(M)) * D;

% Exract submatrices
M1 = zeros(npts-1,npts-1);
M1(1,:) = MD(npts,1:npts-1);

M2 = MD(1:npts-1,1:npts-1);

M3 = zeros(npts-1,npts-1);
M3(:,1) = MD(1:npts-1,npts);

myweights;
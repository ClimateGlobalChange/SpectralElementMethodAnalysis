function [M1,M2,M3] = se_derivative_matrices(npts)
% se_derivative_matrices - Generate matrices for the spectral element
% first derivative operator with unit element width.
%
% Syntax:  [M1,M2,M3] = se_derivative_matrices(npts)
%
% Inputs:
%    npts - Order of the spectral element method
%
% Outputs:
%    M1 - Matrix of left coefficients (npts-1 x npts-1)
%    M2 - Matrix of center coefficients (npts-1 x npts-1)
%    M3 - Matrix of right coefficients (npts-1 x npts-1)
%
% Example (derivative of a line): 
%    >> [M1,M2,M3] = se_derivative_matrices(4);
%    >> [pts,~] = lglnodes(4, 0, 1);
%    >> M1*(pts(1:3)-1)+M2*pts(1:3)+M3*(pts(1:3)+1)
%
%    ans =
%       1.0000
%       1.0000
%       1.0000
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 08-Oct-2017

%------------- BEGIN CODE --------------

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[myroots,~] = lglnodes(npts,0,1);

% Find coefficients of characteristic polynomials
[~,difflegpoly] = characteristic_polyfit(myroots);

% Compute evolution matrices
M1 = zeros(npts-1);
M2 = zeros(npts-1);
M3 = zeros(npts-1);

for m = 1:npts-1
    M1(1,m) = polyval(difflegpoly(m,1:npts), myroots(npts));
end
M1(1,:) = 0.5 * M1(1,:);

for n = 1:npts-1
    for m = 1:npts-1
        M2(m,n) = polyval(difflegpoly(n,:), myroots(m));
    end
end
M2(1,:) = 0.5 * M2(1,:);
M2(1,1) = 0.5 * polyval(difflegpoly(1,1:npts), myroots(1)) ...
        + 0.5 * polyval(difflegpoly(npts,1:npts), myroots(npts));

for n = 1:npts-1
    M3(n,1) = polyval(difflegpoly(npts,1:npts), myroots(n));
end
M3(1,1) = 0.5 * M3(1,1);

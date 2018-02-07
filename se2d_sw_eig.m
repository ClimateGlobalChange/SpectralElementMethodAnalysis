function pevals = se2d_sw_eig(npts, thetax, thetay, M)
% se2d_sw_eig - Compute the eigenvalues of the semi-discrete evolution 
% matrix for the 2D non-dimensionalized linearized shallow-water equations
% with spectral element discretization.  That is, the matrix M in equation
% dphi_j/dt = M phi_j.  Sort the resulting eigenvalues by wavenumber using
% se2d_sw_isolate_physical_mode().
%
% Syntax:  pevals = se2d_sw_eig(npts, thetax, thetay, M)
%
% Inputs:
%    npts - Order of the spectral element method
%    thetax - Values of the non-dimensional wavenumber in the x direction
%             (specified as a vector)
%    thetax - Values of the non-dimensional wavenumber in the y direction
%             (specified as a vector)
%    M - Evolution matrix, for instance from se2d_sw_matrix_exact()
%
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 06-Feb-2018
%
%------------- BEGIN CODE --------------
if (~isscalar(npts))
    error('npts argument must be a scalar');
end
if (~isvector(thetax))
    error('thetax argument must be a vector');
end
if (~isvector(thetay))
    error('thetay argument must be a vector');
end

npx = 3*(npts-1)*(npts-1);

if ((size(M,1) ~= npx) || (size(M,2) ~= npx))
    error('M argument must be a squre matrix of size 3*(npts-1)^2');
end

% Compute all eigenvalues and eigenvectors
evals = zeros(npx,length(thetax));
evecs = zeros(npx,npx,length(thetax));

% Calculate eigenvalues
for j = 1:length(thetax)
    [V,D] = eig(M(:,:,j));
    eigX = diag(D);
    
    evals(:,j) = eigX;
    evecs(:,:,j) = V;
end

% Isolate the physical mode
pevals = se2d_sw_isolate_physical_mode(npts, thetax, thetay, evals, evecs);

end
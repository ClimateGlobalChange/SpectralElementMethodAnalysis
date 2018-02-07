function [evals, evecs] = se_eig(npts, theta, M)
% se_eig - Calculate the eigenvalues, eigenvectors, and identify the
% "physical mode" from the array of semi-discrete evolution matrices.
%
% Syntax:  [evals, evecs] = se_eig(npts, theta, M)
%
% Inputs:
%    npts - Order of the spectral element method
%    theta - Vector of non-dimensional wavenumbers (k \Delta x_e)
%    M - Semi-discrete evolution matrix (npts-1 x npts-1 x length(theta))
%
% Outputs:
%    evals - Eigenvalues of M (npts-1 x length(theta))
%    evecs - Eigenvectors of M (npts-1 x npts-1 x length(theta))
%
%    The eigenvalues and eigenvectors are sorted so that the positive
%    physical mode is first in the array, ordered along with wavenumber.
%
% Example: 
%    >> M = se_advection_matrix_exact(4, 0.1, 4, 0.001);
%    >> [evals, evecs, ix_phys] = se_eig(4, 0.1, M);
%
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 08-Oct-2017

%------------- BEGIN CODE --------------

% Check arguments
if (~isscalar(npts))
    error('npts argument must be a scalar');
end
if (~isvector(theta))
    error('theta argument must be a vector');
end

% Compute all eigenvalues and eigenvectors
evals = zeros(npts-1,length(theta));
evecs = zeros(npts-1,npts-1,length(theta));

for j = 1:length(theta)
    [vec,val] = eig(M(:,:,j));
    evals(:,j) = diag(val);
    evecs(:,:,j) = vec;
end

% Isolate the physical mode
[evals, evecs] = se_isolate_physical_mode(npts, theta, evals, evecs);

end
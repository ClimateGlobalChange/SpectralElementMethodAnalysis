function M = se_advection_matrix_exact(npts, theta, nuorder, nuvalue, nutype)
% se_advection_matrix_exact - Calculate the semi-discrete evolution matrix
% for the 1D non-dimensionalized advection equation with spectral element
% discretization.  That is, the matrix M in equation dphi_j/dt = M phi_j.
%
% Syntax:  M = se_advection_matrix_exact(npts, theta, nuorder, nuvalue, nutype)
%
% Inputs:
%    npts - Order of the spectral element method
%    theta - Vector of non-dimensional wavenumbers (k \Delta x_e)
%    nuorder - Order of diffusion to apply (must be even)
%    nuvalue - Non-dimensional diffusion coefficient
%    nutype - 1 [default] (use laplacian operator for diffusion)
%           - 2 (use derivative operator for diffusion)
%
% Outputs:
%    M - Semi-discrete evolution matrix (npts-1 x npts-1 x length(theta))
%
% Example: 
%    >> M = se_advection_matrix_exact(4, 0.1, 4, 0.001)
%
%    M =
%      -3.0020 + 0.0998i   7.0857 - 0.0952i  -4.0802 + 0.1950i
%      -1.6317 - 0.0853i  -0.7505 - 0.0000i   2.3878 + 0.0335i
%       2.8239 + 0.2451i  -2.0843 - 0.0335i  -0.7505 - 0.0000i
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
if (~exist('nuvalue'))
    if (exist('nuorder'))
        error('Diffusion requires both nuorder and nuvalue to be specified');
    end
    nuorder = 0;
    nuvalue = 0;
end
if (mod(nuorder, 2) ~= 0)
    error('nuorder argument must be even');
end
if (nuvalue < 0)
    disp('WARNING: anti-diffusion being applied in calculation');
end
if (~exist('nutype'))
    nutype = 1;
end
if ((nutype ~= 1) && (nutype ~= 2))
    error('nutype argument must take value 1 or 2');
end

% Generate SE derivative matrices
[M1,M2,M3] = se_derivative_matrices(npts);

% Compute advection operator
M = zeros(npts-1, npts-1, length(theta));

% Without diffusion
if (nuorder == 0)
    for j = 1:length(theta)
    	M(:,:,j) = M1 * exp(-1i * theta(j)) + M2 + M3 * exp(1i * theta(j));
    end

% With laplacian-operator diffusion
elseif (nutype == 1)
    [D1,D2,D3] = se_laplacian_matrices(npts);
    for j = 1:length(theta)
        D = D1 * exp(-1i * theta(j)) + D2 + D3 * exp(1i * theta(j));
    	M(:,:,j) = M1 * exp(-1i * theta(j)) + M2 + M3 * exp(1i * theta(j));
        M(:,:,j) = M(:,:,j) - (-1)^(nuorder/2) * nuvalue * D^(nuorder/2);
    end
    
% With derivative-operator diffusion
elseif (nutype == 2)
    for j = 1:length(theta)
    	M(:,:,j) = M1 * exp(-1i * theta(j)) + M2 + M3 * exp(1i * theta(j));
        M(:,:,j) = M(:,:,j) - (-1)^(nuorder/2) * nuvalue * M(:,:,j)^(nuorder);
    end
end

end

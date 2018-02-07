function MDall = se2d_sw_matrices(npts, As, xonly)
% se2d_sw_matrices - Generate matrices for the spectral element
% 2D non-dimensional linearized shallow water operator.
%
% Syntax:  MDall = se2d_scalar_laplacian_matrices(npts, xonly)
%
% Inputs:
%    npts - Order of the spectral element method
%    As - Ratio of the grid spacing to Rossby radius of deformation
%    xonly - If 1, only include the second derivative term in the x
%            direction (default 0)
%
% Outputs:
%    MDall - Matrix of coefficients for the linearized shallow water
%          operator for all neighboring spectral elements
%          (3*(npts-1) x 3*(npts-1) x 3 x 3)
%
% Remarks:
%    Degrees of freedom within the 2D spectral element are first given in
%    the x direction, then in the y direction.
%
% Example usage:
%    >> MDall = se2d_sw_matrices(4,0.5,0);
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 06-Feb-2018

%------------- BEGIN CODE --------------
if (~isscalar(npts))
    error('argument npts must be of type scalar');
end
if (~isscalar(As))
    error('argument As must be of type scalar');
end
if (~exist('xonly'))
    xonly = 0;
end
if ((xonly ~= 0) && (xonly ~= 1))
    error('invalid value for agument xonly; value should be 0 or 1');
end

% Get the 2D first derivative matrices in the x and y directions
[MDx, MDy] = se2d_deriv_matrices(npts);

% Spectral element size
npx = (npts-1)*(npts-1);

% Extract all submatrices
MDall = zeros(3*npx, 3*npx, 3, 3);

% Height update equation
MDall(1:npx, npx+1:2*npx, :, :) = MDall(1:npx, npx+1:2*npx, :, :) - MDx;
if (xonly ~= 1)
    MDall(1:npx, 2*npx+1:3*npx, :, :) = MDall(1:npx, 2*npx+1:3*npx, :, :) - MDy;
end

% X velocity update equation
MDall(npx+1:2*npx, 1:npx, :, :) = MDall(npx+1:2*npx, 1:npx, :, :) - MDx;
MDall(npx+1:2*npx, 2*npx+1:3*npx, 2, 2) = MDall(npx+1:2*npx, 2*npx+1:3*npx, 2, 2) + As*eye(npx);

% Y velocity update equation
if (xonly ~= 1)
    MDall(2*npx+1:3*npx, 1:npx, :, :) = MDall(2*npx+1:3*npx, 1:npx, :, :) - MDy;
end
MDall(2*npx+1:3*npx, npx+1:2*npx, 2, 2) = MDall(2*npx+1:3*npx, npx+1:2*npx, 2, 2) - As*eye(npx);


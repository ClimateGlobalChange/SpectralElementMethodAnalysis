function Mevol = se2d_sw_matrix_exact(npts, theta, thetadir, As, nuorder, InvReS, InvReD, InvReV, xonly)
% se2d_sw_matrix_exact - Generate the matrix for the spectral element
% non-dimensional linearized shallow water operator with hyperdiffusion.
%
% Syntax:  Mevol = se2d_sw_matrix_exact(npts, theta, thetadir, As, nuorder, InvReS, InvReD, InvReV, xonly)
%
% Inputs:
%    npts - Order of the spectral element method
%    theta - Values of the non-dimensional wavenumber in the x direction
%    thetadir - Angle of wave propagation (used to calculate thetay)
%    As - Ratio of the grid spacing to Rossby radius of deformation
%    nuorder - Order of diffusion to apply (must be even)
%    InvReS - Inverse scalar hyper-Reynolds number
%    InvReD - Inverse divergence damping hyper-Reynolds number
%    InvReV - Inverse vorticity damping hyper-Reynolds number
%    xonly - Ignore all derivative terms in the y direction
%
% Remarks:
%    Degrees of freedom within the 2D spectral element are first given in
%    the x direction, then in the y direction. 
%
% Example usage:
%    >> Mevol = se2d_sw_matrix_exact(4, [0:0.1:1]*pi, 0, 0.5, 4, 1e-4, 1e-4, 1e-4, 0);
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 06-Feb-2018
%
%------------- BEGIN CODE --------------

if (~isscalar(npts))
    error('argument npts must be of type scalar');
end
if (~isvector(theta))
    error('argument theta must be of type vector');
end
if ((thetadir < -pi/4) || (thetadir > pi/4))
    error('thetadir out of range [-pi/4,pi/4]');
end
if (~exist('xonly'))
    xonly = 0;
end
if ((xonly ~= 0) && (xonly ~= 1))
    error('invalid value for agument xonly; value should be 0 or 1');
end
if ((xonly == 1) && (thetadir ~= 0))
    error('if xonly is set, thetadir must be 0');
end

% Get the shallow water update matrices
MDall = se2d_sw_matrices(npts, As, xonly);

% Get the 2D scalar hyperdiffusion matrices
MDdiff = se2d_scalar_laplacian_matrices(npts, xonly);

% Get the 2D vector hyperdiffusion matrices
if (xonly == 0)
    [DDall_uu, DDall_uv, DDall_vu, DDall_vv, VDall_uu, VDall_uv, VDall_vu, VDall_vv] = se2d_vecdiff_matrices(npts);
end

% Spectral element size
npx = (npts-1)*(npts-1);

% Generate thetas in each direction
thetax = theta;
thetay = theta * tan(thetadir);

% Transformed matrix
Mevol = zeros(3*npx, 3*npx, length(theta));

for j = 1:length(theta)
    ScaDamp = zeros(npx,npx);
    for n = 1:3
    for m = 1:3
        Mevol(:,:,j) = Mevol(:,:,j) + MDall(:,:,n,m) * exp(1i * (m-2) * thetax(j)) * exp(1i * (n-2) * thetay(j));
        ScaDamp = ScaDamp + MDdiff(:,:,n,m) * exp(1i * (m-2) * thetax(j)) * exp(1i * (n-2) * thetay(j));
    end
    end

    % Add scalar diffusion to height field
    Mevol(1:npx,1:npx,j) = Mevol(1:npx,1:npx,j) - (-1)^(nuorder/2) * InvReS * ScaDamp^(nuorder/2);

    % Add divergence damping to velocity field
    if (xonly == 1)
        Mevol(2*npx+1:3*npx,2*npx+1:3*npx,j) = Mevol(1:npx,1:npx,j) - (-1)^(nuorder/2) * InvReD * ScaDamp^(nuorder/2);

    else
        % Divergence and vorticity damping
        DivDamp = zeros(3*npx, 3*npx);
        VorDamp = zeros(3*npx, 3*npx);
        for n = 1:3
        for m = 1:3
            waveform = exp(1i * (m-2) * thetax(j)) * exp(1i * (n-2) * thetay(j));

            DivDamp(npx+1:2*npx,npx+1:2*npx) = DivDamp(npx+1:2*npx,npx+1:2*npx) + DDall_uu(:,:,n,m) * waveform;
            DivDamp(npx+1:2*npx,2*npx+1:3*npx) = DivDamp(npx+1:2*npx,2*npx+1:3*npx) + DDall_uv(:,:,n,m) * waveform;
            DivDamp(2*npx+1:3*npx,npx+1:2*npx) = DivDamp(2*npx+1:3*npx,npx+1:2*npx) + DDall_vu(:,:,n,m) * waveform;
            DivDamp(2*npx+1:3*npx,2*npx+1:3*npx) = DivDamp(2*npx+1:3*npx,2*npx+1:3*npx) + DDall_vv(:,:,n,m) * waveform;

            VorDamp(npx+1:2*npx,npx+1:2*npx) = VorDamp(npx+1:2*npx,npx+1:2*npx) + VDall_uu(:,:,n,m) * waveform;
            VorDamp(npx+1:2*npx,2*npx+1:3*npx) = VorDamp(npx+1:2*npx,2*npx+1:3*npx) + VDall_uv(:,:,n,m) * waveform;
            VorDamp(2*npx+1:3*npx,npx+1:2*npx) = VorDamp(2*npx+1:3*npx,npx+1:2*npx) + VDall_vu(:,:,n,m) * waveform;
            VorDamp(2*npx+1:3*npx,2*npx+1:3*npx) = VorDamp(2*npx+1:3*npx,2*npx+1:3*npx) + VDall_vv(:,:,n,m) * waveform;
        end
        end

        Mevol(:,:,j) = Mevol(:,:,j) - (-1)^(nuorder/2) * InvReD * DivDamp^(nuorder/2);
        Mevol(:,:,j) = Mevol(:,:,j) - (-1)^(nuorder/2) * InvReV * VorDamp^(nuorder/2);
    end
end



%MDall(1:npx, 1:npx, :, :) = MDall(1:npx, 1:npx, :, :) + InvReS * MDdiff;
%MDall(npx+1:2*npx, npx+1:2*npx, :, :) = MDall(npx+1:2*npx, npx+1:2*npx, :, :) + InvReV * MDdiff;
%MDall(2*npx+1:3*npx, 2*npx+1:3*npx, :, :) = MDall(2*npx+1:3*npx, 2*npx+1:3*npx, :, :) + InvReV * MDdiff;

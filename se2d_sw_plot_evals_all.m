function se2d_sw_plot_evals_all(npts, thetadir, As, nuorder, InvReS, InvReD, InvReV, xonly)
% se2d_sw_plot_evals_all - Plot the eigenvalues of the semi-discrete 
% evolution matrix for the 2D non-dimensionalized linearized shallow-water
% equations with spectral element discretization.  That is, the matrix M in equation
% dphi_j/dt = M phi_j.  Sorting by wavenumber is not performed.
%
% Syntax:  se2d_sw_plot_evals(npts, thetadir, As, nuorder, InvReS, InvReD, InvReV, xonly)
%
% Inputs:
%    npts - Order of the spectral element method
%    thetadir - Direction of wave propagation (-pi/4 <= thetadir <= pi/4)
%    As - Ratio of the grid spacing to Rossby radius of deformation
%    nuorder - Order of diffusion to apply (must be even)
%    InvReS - Inverse scalar hyper-Reynolds number
%    InvReD - Inverse divergence damping hyper-Reynolds number
%    InvReV - Inverse vorticity damping hyper-Reynolds number
%    xonly - Neglect all derivative terms in the y direction
%
% Example: 
%    >> se2d_sw_plot_evals_all(4, pi/4, 0.5, 4, 1e-4, 1e-4, 1e-4, 0);
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
if (~exist('As'))
    As = 0;
end
if (~isscalar(As))
    error('As argument must be a scalar');
end
if (~exist('nuorder'))
    nuorder = 0;
    InvReS = 0;
    InvReD = 0;
    InvReV = 0;
end
if (mod(nuorder, 2) ~= 0)
    error('nuorder argument must be even');
end
if ((~exist('InvReS')) || (~exist('InvReD')) || (~exist('InvReV')))
    if (exist('nuorder'))
        error('Diffusion requires both nuorder and InvReS, InvReD, InvReV to be specified');
    end
    nuorder = 0;
    InvReS = 0;
    InvReD = 0;
    InvReV = 0;
end
if (~exist('xonly'))
    xonly = 0;
end

% Wave modes to plot
theta = [0:0.001:1] * 2 * pi;

% Generate thetas in each direction
thetax = theta;
thetay = theta * tan(thetadir);

% Get the evolution matrices
Mevol = se2d_sw_matrix_exact(npts, theta, thetadir, As, nuorder, InvReS, InvReD, InvReV, xonly);

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[gll,gwt] = lglnodes(npts,0,1);
gwt2d = gwt(1:npts-1) * gwt(1:npts-1)';
W = diag(reshape(gwt2d,[1 (npts-1)*(npts-1)]));

% Calculate eigenvalues
alleigs = zeros(length(theta),3*(npts-1)*(npts-1));
ieigs = zeros(length(theta),3*(npts-1)*(npts-1));
reigs = zeros(length(theta),3*(npts-1)*(npts-1));
for j = 1:length(theta)
    [V,D] = eig(Mevol(:,:,j));
    eigX = diag(D);

    for k = 1:3*(npts-1)*(npts-1)
        alleigs(j,k) = eigX(k);
    end

    alleigs(j,:) = (sortrows([imag(alleigs(j,:).') real(alleigs(j,:).')])*[1i;1]).';
    ieigs(j,:) = sort(imag(alleigs(j,:).'));
    reigs(j,:) = sort(real(alleigs(j,:).'));
end

% Plot
disp(sprintf('Stability criteria (should be 0): %1.5e', max(max(reigs))));

width = 2;

gcf = figure(1);
set(gcf, 'Position', [100 100 1200 600]); 

axis1 = axes('Position', [0.08 0.2 0.4 0.65]);
plot(theta, ieigs, '-', 'LineWidth', 2);
minr = min(min(ieigs));
maxr = max(max(ieigs));
if (minr == maxr)
    axis([0 width*pi -0.1 0.1]);
else
    axis([0 width*pi minr-0.05*(maxr-minr) maxr+0.05*(maxr-minr)]);
end
set(axis1, 'FontSize', 16);
xlabel('Dimensionless wavenumber (k_e \Delta x_e)');
ylabel('\omega_i');

axis2 = axes('Position', [0.58 0.2 0.4 0.65]);
plot(theta, reigs, '-', 'LineWidth', 2);
minr = min(min(reigs));
maxr = max(max(reigs));
if (minr == maxr)
    axis([0 width*pi -0.1 0.1]);
else
    axis([0 width*pi minr-0.05*(maxr-minr) maxr+0.05*(maxr-minr)]);
end
set(axis2, 'FontSize', 16);
xlabel('Dimensionless wavenumber (k_e \Delta x_e)');
ylabel('\omega_r');

%hold on;
%plot(theta, sqrt(InvRo^2 + theta.^2 + (theta .* tan(thetadir)).^2), 'k-', 'LineWidth', 1);
%plot(theta, sqrt(InvRo^2 + theta.^2 + (2*pi)^2 + (theta .* tan(thetadir)).^2), 'k-', 'LineWidth', 1);
%plot(theta, 0*theta, 'k-', 'LineWidth', 1);
%plot(theta, -sqrt(InvRo^2 + theta.^2 + (theta .* tan(thetadir)).^2), 'k-', 'LineWidth', 2);
%hold off;

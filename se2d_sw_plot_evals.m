function [theta_ext, alleigs_out] = se2d_sw_plot_evals(npts, As, nuorder, InvReS, InvReD, InvReV)
% se2d_sw_plot_evals - Plot the eigenvalues of the semi-discrete evolution 
% matrix for the 2D non-dimensionalized linearized shallow-water equations
% with spectral element discretization.  That is, the matrix M in equation
% dphi_j/dt = M phi_j.  Eigenvalues are sorted by wavenumber using 
% se2d_sw_isolate_physical_mode().
%
% Syntax:  se2d_sw_plot_evals(npts, As, nuorder, InvReS, InvReD, InvReV)
%
% Inputs:
%    npts - Order of the spectral element method
%    As - Ratio of the grid spacing to Rossby radius of deformation
%    nuorder - Order of diffusion to apply (must be even)
%    InvReS - Inverse scalar hyper-Reynolds number
%    InvReD - Inverse divergence damping hyper-Reynolds number
%    InvReV - Inverse vorticity damping hyper-Reynolds number
%
% Example: 
%    >> se2d_sw_plot_evals(4, 0.5, 4, 1e-4, 1e-4, 0);
%
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 02-Feb-2018
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

% Set thetadir to zero
thetadir = 0;

% xonly not implemented
if (exist('xonly'))
    error('xonly not implemented');
end
xonly = 0;

% Wave modes to plot
theta = [0:0.001:1] * 2 * pi;

% Generate thetas in each direction
thetax = theta;
thetay = theta * tan(thetadir);

% Get the evolution matrices
Mevol = se2d_sw_matrix_exact(npts, theta, thetadir, As, nuorder, InvReS, InvReD, InvReV, xonly);

% Calculate physical eigenvalues
pevals = se2d_sw_eig(npts, thetax, thetay, Mevol);

% Extend theta
theta_ext = zeros([1 length(theta)*(npts-1)]);
for j = 1:length(theta)
    for m = 1:npts-1
        theta_ext((m-1)*length(theta)+j) = theta(j) + 2*pi*(m-1);
    end
end

% Range of frequencies
mini = min(min(min(imag(pevals))));
maxi = max(max(max(imag(pevals))));

minr = min(min(min(real(pevals))));
maxr = max(max(max(real(pevals))));

% Output stability criteria
fprintf('Stability criteria (should be <= 0): %1.5e\n', maxr);

% Current figure
gcf = figure(1);

% Loop through all modes in the y direction
for t = 1:(npts-1)
    alleigs = reshape(pevals(:,t,:), [3 (npts-1)*length(theta)]);

    for j = 1:(npts-1)*length(theta)
        alleigs(:,j) = sortrows([imag(alleigs(:,j)) real(alleigs(:,j))])*[1i;1];
    end
    reigs = real(alleigs);
    ieigs = imag(alleigs);

    if (t == 1)
        alleigs_out = alleigs;
        
        set(gcf, 'Position', [100 100 1200 600]); 

        width = (npts-1);

        axis1 = axes('Position', [0.08 0.2 0.4 0.65]);
        plot(theta_ext, ieigs, 'k.', 'LineWidth', 2);
        if (mini == maxi)
            axis(axis1, [0 width*pi -0.1 0.1]);
        else
            axis(axis1, [0 width*pi mini-0.05*(maxi-mini) maxi+0.05*(maxi-mini)]);
        end
        set(axis1, 'FontSize', 16);
        title('Imaginary [\lambda] (Frequency)');
        xlabel('Dimensionless wavenumber (k_e \Delta x_e)');
        ylabel('\omega_i');

        axis2 = axes('Position', [0.58 0.2 0.4 0.65]);
        plot(theta_ext, reigs, 'k.', 'LineWidth', 2);
        if (minr == maxr)
            axis(axis2, [0 width*pi -0.1 0.1]);
        else
            axis(axis2, [0 width*pi minr-0.05*(maxr-minr) maxr+0.05*(maxr-minr)]);
        end
        set(axis2, 'FontSize', 16);
        title('Real [\lambda] (Diffusivity)');
        xlabel('Dimensionless wavenumber (k_e \Delta x_e)');
        ylabel('\omega_r');
    else
        fmt = 'r.';
        if (t == 3)
            fmt = 'g.';
        elseif (t == 4)
            fmt = 'm.';
        elseif (t == 5)
            fmt = 'b.';
        elseif (t == 6)
            fmt = 'c.';
        end
        hold(axis1, 'on');
        plot(axis1, theta_ext, ieigs, fmt, 'LineWidth', 2);
        hold(axis1, 'off');
        hold(axis2, 'on');
        plot(axis2, theta_ext, reigs, fmt, 'LineWidth', 2);
        hold(axis2, 'off');
    end

end


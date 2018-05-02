function [theta_ext, pevals] = se_plot_evals(npts, nuorder, nuvalue, nutype, theta_numpts)
% se_plot_evals - Plot the eigenvalues of the semi-discrete evolution matrix
% for the 1D non-dimensionalized advection equation with spectral element
% discretization.  That is, the matrix M in equation dphi_j/dt = M phi_j.
%
% Syntax:  se_plot_evals(npts, nuorder, nuvalue, nutype, theta_numpts)
%
% Inputs:
%    npts - Order of the spectral element method
%    nuorder - Order of diffusion to apply (must be even)
%    nuvalue - Non-dimensional diffusion coefficient
%    nutype - 1 [default] (use norder/2 applications of laplacian operator for diffusion)
%           - 2 (use norder applications of derivative operator for diffusion)
%    theta_numpts - Number of values of theta to use in the plot (default 10000)
%
% Example: 
%    >> se_plot_evals(4, 4, 0.001);
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

if (~exist('theta_numpts'))
    theta_numpts = 1000;
else
    theta_numpts = floor(theta_numpts);
end
if (~isscalar(theta_numpts))
    error('theta_numpts argument must be a scalar');
end
if (theta_numpts < 2)
    error('theta_numpts argument must be at least 2');
end

% uniformly distribute theta
theta = [0:theta_numpts-1]/(theta_numpts-1) * 2 * pi;

% Advection matrices
M = se_advection_matrix_exact(npts, theta, nuorder, nuvalue, nutype);

% Get evolution matrices
[evals, ~] = se_eig(npts, theta, M);

% Extend theta
theta_ext = zeros([1 length(theta)*(npts-1)]);
for j = 1:length(theta)
    for m = 1:npts-1
        theta_ext((m-1)*length(theta)+j) = theta(j) + 2*pi*(m-1);
    end
end
pevals = reshape(evals.', [1 length(theta)*(npts-1)]);

% Break into imaginary and real components
ievals = imag(evals.');
revals = real(evals.');

for j = 1:length(theta)
    ievals(j,:) = sort(ievals(j,:));
    revals(j,:) = sort(revals(j,:));
end

evals = (revals + 1i * ievals).';

% Plot type
width = (npts-1);

% Plot frequency (imaginary component of eigenvalues)
axis1 = axes('Position', [0.08 0.2 0.4 0.65]);
%plot(theta, ievals, 'k--', 'LineWidth', 1, 'Parent', axis1);
plot([0 width*pi], [0 width*pi], 'k-');
hold on;
plot(theta_ext, imag(pevals), 'k-', 'LineWidth', 2);
hold off;
set(gca, 'FontSize', 14);
xlabel('Non-dim. wavenumber (\theta)');
title('Frequency (\omega_i)');

plotwidth = max(max(ievals)) - min(min(ievals));
axis([0 width*pi min(min(ievals))-0.05*plotwidth max(max(ievals))+0.05*plotwidth]);

% Plot diffusion (real component of eigenvalues)
axis2 = axes('Position', [0.58 0.2 0.4 0.65]);
%plot(theta, revals, 'k--', 'Parent', axis2);
%hold on;

plot(theta_ext, real(pevals), 'k-', 'LineWidth', 2);
%hold off;
set(gca, 'FontSize', 14);
xlabel('Non-dim. wavenumber (\theta)');
title('Diffusivity (\omega_r)');

plotwidth = max(max(ievals)) - min(min(ievals));
axis([0 width*pi min(min(revals))-0.05*plotwidth 0.1]);

% Add title
antitle = annotation('textbox', [0.08 0.9 0.84 0.08]);
set(antitle, 'String', 'SE eigenstructure');
set(antitle, 'FontSize', 16);
set(antitle, 'HorizontalAlignment', 'center');
set(antitle, 'LineStyle', 'none');

end

function se_plot_evals_dt(npts, temporal_scheme, cfl, nuorder, nuvalue, nutype, numethod, theta_numpts)
% se_plot_evals - Plot the eigenvalues of the fully discrete evolution matrix
% for the 1D non-dimensionalized advection equation with spectral element
% discretization.
%
% Syntax:  se_plot_evals_dt(npts, temporal_scheme, cfl, nuorder, nuvalue, nutype)
%
% Inputs:
%    npts - Order of the spectral element method
%    temporal_scheme - Temporal integration scheme from time_discrete_M
%    cfl - Non-dimensional Courant number to use (c * dt / dx), where
%          dx is the average distance between degrees of freedom
%    nuorder - Order of diffusion to apply (must be even)
%    nuvalue - Non-dimensional diffusion coefficient
%    nutype - 1 (use norder applications of derivative operator for diffusion)
%             2 [default] (use norder/2 applications of laplacian operator for diffusion)
%    numethod - 'inline' (apply hyperdiffusion as part of each RK subcycle)
%             - 'split' (apply hyperdiffusion with FE method after all RK subcycles)
%    theta_numpts - Number of values of theta to use in the plot (default 10000)
%
% Example: 
%    >> se_plot_evals_dt(4, 'KG53', 0.5, 4, 0.001, 1, 'split');
%
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 08-Oct-2017

%------------- BEGIN CODE --------------


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
if (~exist('numethod'))
    numethod = 'inline';
end
if (strcmp(numethod, 'inline') && strcmp(numethod, 'split'))
    error('numethod must be one of inline or split');
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

% Rescale the Courant number to use average spacing between degrees of
% freedom.
cfl = cfl / (npts-1);

% uniformly distribute theta
theta = [0:theta_numpts-1]/(theta_numpts-1) * 2 * pi;

% Generate SE derivative matrices
[D1,D2,D3] = se_derivative_matrices(npts);

% Compute advection and hyperdiffusion operators
Dx = zeros(npts-1, npts-1, length(theta));
HDiff = zeros(npts-1, npts-1, length(theta));

% With derivative-operator diffusion
if (nutype == 2)
    for j = 1:length(theta)
    	Dx(:,:,j) = D1 * exp(-1i * theta(j)) + D2 + D3 * exp(1i * theta(j));
        HDiff(:,:,j) = - (-1)^(nuorder/2) * nuvalue * Dx(:,:,j)^(nuorder);
    end

% With laplacian-operator diffusion
elseif (nutype == 1)
    [L1,L2,L3] = se_laplacian_matrices(npts);
    for j = 1:length(theta)
    	Dx(:,:,j) = D1 * exp(-1i * theta(j)) + D2 + D3 * exp(1i * theta(j));
        HDiff(:,:,j) = L1 * exp(-1i * theta(j)) + L2 + L3 * exp(1i * theta(j));
        HDiff(:,:,j) = - (-1)^(nuorder/2) * nuvalue * HDiff(:,:,j)^(nuorder/2);
    end
end

% Compute discrete advection update operator
M = zeros(npts-1, npts-1, length(theta));

% Inline diffusion
if (strcmp(numethod, 'inline'))
    for j = 1:length(theta)
        M(:,:,j) = time_discrete_M(temporal_scheme, Dx(:,:,j) + HDiff(:,:,j), cfl);
    end
end
if (strcmp(numethod, 'split'))
    for j = 1:length(theta)
        M(:,:,j) = time_discrete_M(temporal_scheme, Dx(:,:,j), cfl);
        M(:,:,j) = (eye(npts-1) + cfl * HDiff(:,:,j)) * M(:,:,j);
    end
end

% Compute eigenvalues
evals = zeros(npts-1, length(theta));
evecs = zeros(npts-1, npts-1, length(theta));
for j = 1:length(theta)
    [evals(:,j),evecs(:,:,j)] = time_discrete_eig(cfl, M(:,:,j));
end

% Isolate the physical mode
[evals, evecs] = se_isolate_physical_mode(npts, theta, evals, evecs);

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
set(gca, 'FontSize', 16);
xlabel('Dimensionless wavenumber (k_e \Delta x_e)');
ylabel('\omega_r');

plotwidth = max(max(ievals)) - min(min(ievals));
axis([0 width*pi min(min(ievals))-0.05*plotwidth max(max(ievals))+0.05*plotwidth]);

% Plot diffusion (real component of eigenvalues)
axis2 = axes('Position', [0.58 0.2 0.4 0.65]);
%plot(theta, revals, 'k--', 'Parent', axis2);
%hold on;

plot(theta_ext, real(pevals), 'k-', 'LineWidth', 2);
%hold off;
set(gca, 'FontSize', 16);
xlabel('Dimensionless wavenumber (k_e \Delta x_e)');
ylabel('\omega_i');

plotwidth = max(max(ievals)) - min(min(ievals));
axis([0 width*pi min(min(revals))-0.05*plotwidth 0.1]);

% Add title
antitle = annotation('textbox', [0.08 0.9 0.84 0.08]);
set(antitle, 'String', 'SE eigenstructure');
set(antitle, 'FontSize', 16);
set(antitle, 'HorizontalAlignment', 'center');
set(antitle, 'LineStyle', 'none');

end

% Convert the time-discrete eigenvalues to values that are more consistent
% with what would be expected from the eigenvalues of the standalone spatial
% discretization.
function [evals,evecs] = time_discrete_eig(cfl, T)

    % Compute eigenvalues
    [evecs,ev] = eig(-T);
    evals = diag(ev);

    for k = 1:size(evals,1)
        evals(k) = - 1i * 1/cfl * atan2(imag(evals(k)), -real(evals(k))) + log(abs(evals(k))^(1/cfl));
    end
end
function [theta,waveerr,cerr,differr] = se_wave_accuracy(npts, nuorder, nuvalue, nutype)
% se_wave_accuracy - Plot wave error, phase speed error, and diffusive
% error for the 1D non-dimensionalized advection equation with spectral
% element discretization.
%
% Syntax:  se_wave_accuracy(npts, nuorder, nuvalue, nutype)
%
% Inputs:
%    npts - Order of the spectral element method
%    nuorder - Order of diffusion to apply (must be even)
%    nuvalue - Non-dimensional diffusion coefficient
%    nutype - 1 [default] (use norder/2 applications of laplacian operator for diffusion)
%           - 2 (use norder applications of derivative operator for diffusion)
%
% Example: 
%    >> se_wave_accuracy(4, 4, 1e-4, 1);
%
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 1-May-2018

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

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[gll,gllw] = lglnodes(npts,0,1);
W = diag(gllw(1:npts-1));

% Number of points
theta_numpts = 100;

% uniformly distribute theta
theta = [0:theta_numpts-1]/(theta_numpts-1) * 2 * pi;

% Advection matrices
M = se_advection_matrix_exact(npts, theta, nuorder, nuvalue, nutype);

% Get evolution matrices
[evals, evecs] = se_eig(npts, theta, M);

% Extend theta
theta_ext = zeros([1 length(theta)*(npts-1)]);
for j = 1:length(theta)
    for m = 1:npts-1
        theta_ext((m-1)*length(theta)+j) = theta(j) + 2*pi*(m-1);
    end
end
pevals = reshape(evals.', [1 length(theta)*(npts-1)]);
pevecs = reshape(permute(evecs, [1 3 2]), [(npts-1) length(theta)*(npts-1)]);

% Compute wave error, phase speed error, and diffusion error
waveerr = 100*ones(size(theta_ext));
cerr = 100*ones(size(theta_ext));
differr = 100*ones(size(theta_ext));

for j = 1:length(theta_ext)
    exactev = exp(1i * theta_ext(j) * gll(1:npts-1));
    discrev = pevecs(:,j).';

    exactev_mag = sqrt(abs(exactev'*W*exactev));
    discrev_mag = sqrt(abs(discrev*W*discrev'));
    if (discrev_mag == 0.0)
        error('Invalid discrete eigenvector');
    end
 
    waveerr(j) = 1 - abs(discrev*W*conj(exactev)) / (discrev_mag * exactev_mag);
    cerr(j) = 1.0 - imag(pevals(j))/theta_ext(j);
    differr(j) = 1.0 - exp(real(pevals(j)));
end

cerr(1) = 0.0;

figure(1);

% Log axes
logabswaveerr = log10(abs(waveerr));
logabscerr = log10(abs(cerr));
logabsdifferr = log10(abs(differr));

% Plot type
width = (npts-1);

% Wave error
axis1 = axes('Position', [0.08 0.2 0.26 0.65]);
plot(theta_ext, logabswaveerr, 'k-', 'LineWidth', 2);
set(gca, 'FontSize', 14);
%xlabel('Non-dim. wavenumber (\theta)');
title('log_{10} Wave Error');

plotwidth = max(max(logabswaveerr)) - min(min(logabswaveerr));
plotmax = max(max(logabswaveerr))+0.05*plotwidth;
axis([0 width*pi -15 0]);

% Phase error
axis2 = axes('Position', [0.38 0.2 0.26 0.65]);
plot(theta_ext, logabscerr, 'k-', 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlabel('Non-dim. wavenumber (\theta)');
title('log_{10} Phase Error');

plotwidth = max(max(logabscerr)) - min(min(logabscerr));
plotmax = max(max(logabscerr))+0.05*plotwidth;
axis([0 width*pi -15 0]);

% Diffusive error
axis2 = axes('Position', [0.68 0.2 0.26 0.65]);
plot(theta_ext, logabsdifferr, 'k-', 'LineWidth', 2);
set(gca, 'FontSize', 14);
%xlabel('Non-dim. wavenumber (\theta)');
title('log_{10} Diffusive Error');

axis([0 width*pi -15 0]);

end
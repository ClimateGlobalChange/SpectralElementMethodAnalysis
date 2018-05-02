function max_courant = se_cfl_condition(npts, temporal_scheme, nuorder, nuvalue, nutype, numethod)
% se_cfl_condition - Determine the CFL condition (maximum stable Courant
% number) for the spectral element discretization of the advection
% equation with diffusion.
%
% Syntax:  se_cfl_condition(npts, temporal_scheme, nuorder, nuvalue, nutype, numethod)
%
% Inputs:
%    npts - Order of the spectral element method
%    temporal_scheme - Temporal integration scheme from time_discrete_M
%    nuorder - Order of diffusion to apply (must be even)
%    nuvalue - Non-dimensional diffusion coefficient (scalar or vector)
%    nutype - 1 [default] (use laplacian operator for diffusion)
%           - 2 (use derivative operator for diffusion)
%           - 3 (use Legendre filter applied to shortest wavelength mode)
%    numethod - 'inline' [default] (apply diffusion at each RK substage)
%             - 'split' (apply diffusion after RK cycle)
%
% Outputs:
%    max_courant - Maximum Courant number permitted by the CFL condition.
%
% Example: 
%    >> se_cfl_condition(4, 'RK3', 4, 0.001);
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
if (~isvector(nuvalue))
    error('nuvalue argument must be a vector');
end
if (sum(nuvalue < 0) > 0)
    disp('WARNING: anti-diffusion being applied in calculation');
end
if (~exist('nutype'))
    nutype = 2;
end
if ((nutype ~= 1) && (nutype ~= 2) && (nutype ~= 3))
    error('nutype argument must take value 1, 2, or 3');
end
if (~exist('numethod'))
    numethod = 'inline';
end
if (strcmp(numethod, 'inline') && strcmp(numethod, 'split'))
    error('numethod must be one of inline or split');
end

% Initialize maximum Courant number
max_courant = zeros(size(nuvalue));

% theta is chosen based on search precision
theta = [0:0.001:1] * pi;

% Loop over all nuvalues
for n = 1:length(nuvalue)

    % Bisection search to find stable courant number
    search_range = [0 3];
    while(1)
        % Current CFL number to test
        dt = 0.5 * sum(search_range);
        stable = 1;

        % Inline diffusion
        if (strcmp(numethod, 'inline'))
            if ((nutype == 1) || (nutype == 2))
                M = se_advection_matrix_exact(npts, theta, nuorder, nuvalue(n), nutype);
                for j = 1:length(theta)
                    M(:,:,j) = time_discrete_M(temporal_scheme, M(:,:,j), dt);
                end
                
            elseif (nutype == 3)
                legcoeffs = ones(npts,1);
                legcoeffs(npts) = legcoeffs(npts) - nuvalue;
                [B1,B2,B3] = se_legfilter_matrices(npts, legcoeffs);

                M = se_advection_matrix_exact(npts, theta);
                for j = 1:length(theta)
                    F(:,:,j) = B1 * exp(-1i * theta(j)) + B2 + B3 * exp(1i * theta(j));
                    M(:,:,j) = (F(:,:,j) - eye(npts-1)) / dt + M(:,:,j);
                    M(:,:,j) = time_discrete_M(temporal_scheme, M(:,:,j), dt);
                end
            end
        end

        % Time split
        if (strcmp(numethod, 'split'))
            M = se_advection_matrix_exact(npts, theta);
            for j = 1:length(theta)
                M(:,:,j) = time_discrete_M(temporal_scheme, M(:,:,j), dt);
            end

            % Diffusion via derivative operator
            if (nutype == 1)
                [M1,M2,M3] = se_derivative_matrices(npts);
                for j = 1:length(theta)
                    MM = M1 * exp(-1i * theta(j)) + M2 + M3 * exp(1i * theta(j));
                    M(:,:,j) = (eye(npts-1) - dt * (-1)^(nuorder/2) * nuvalue(n) * MM^(nuorder)) * M(:,:,j);
                end
            end

            % Diffusion via Laplacian operator
            if (nutype == 2)
                [D1,D2,D3] = se_laplacian_matrices(npts);
                for j = 1:length(theta)
                    D = D1 * exp(-1i * theta(j)) + D2 + D3 * exp(1i * theta(j));
                    M(:,:,j) = (eye(npts-1) - dt * (-1)^(nuorder/2) * nuvalue(n) * D^(nuorder/2)) * M(:,:,j);
                end
            end

            % Diffusion via Boyd-Vandeven filter
            if (nutype == 3)
                lfcoeffs = ones(npts,1);
                lfcoeffs(npts) = lfcoeffs(npts) - nuvalue;
                [B1,B2,B3] = se_legfilter_matrices(npts, lfcoeffs);
                for j = 1:length(theta)
                    F(:,:,j) = B1 * exp(-1i * theta(j)) + B2 + B3 * exp(1i * theta(j));
                    M(:,:,j) = F(:,:,j) * M(:,:,j);
                end
            end
        end

        % Check stability across all wavenumbers
        for j = 1:length(theta)

            % Compute eigenvalues
            evals = eig(M(:,:,j));

            for k = 1:length(evals)
                if (abs(evals(k)) > 1.0 + 1.0e-14)
                    stable = 0;
                    break;
                end
            end

            if (stable == 0)
                break;
            end
        end

        if (stable == 1)
            search_range(1) = dt;
        else
            search_range(2) = dt;
        end

        if (search_range(2) - search_range(1) < 0.00001)
            break;
        end
    end

    % Rescale maximum Courant number to the average distance between
    % degrees of freedom.
    max_courant(n) = (npts-1) * 0.5 * sum(search_range);   
end

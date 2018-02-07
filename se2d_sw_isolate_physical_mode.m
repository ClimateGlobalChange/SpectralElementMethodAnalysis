function pevals = se2d_sw_isolate_physical_mode(npts, thetax, thetay, evals, evecs)
% se2d_sw_isolate_physical_mode - Sort the eigenvalues in order of 
% increasing physical wavenumber, effectively mapping the wave modes to
% their physically relevant counterparts.
%
% Syntax:  pevals = se2d_sw_isolate_physical_mode(npts, thetax, thetay, evals, evecs)
%
% Inputs:
%    npts - Order of the spectral element method
%    thetax - Vector of non-dimensional wavenumbers in the x direction
%    thetay - Vector of non-dimensional wavenumbers in the y direction
%    evals - Set of eigenvalues (3*(npts-1)^2 x length(theta))
%    evecs - Set of eigenvectors (3*(npts-1)^2 x 3*(npts-1)^2 x length(theta))
%
% Output:
%    pevals - Sorted eigenvalues (3, (npts-1), (npts-1)*length(thetax)),
%      ordered by field, y-wavenumber, and physical x-direction wavenumber.
%
% Remarks:
%    This function is normally called from se2d_sw_eig.
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 06-Feb-2018

%------------- BEGIN CODE --------------

% Output array
pevals = zeros(3,(npts-1),(npts-1)*length(thetax));
pevecs = zeros(3,(npts-1),(npts-1)*length(thetax),(npts-1));

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[gll,gwt] = lglnodes(npts,0,1);
gwt(1) = gwt(1) * 2;
W = diag(reshape(gwt(1:npts-1)*gwt(1:npts-1)',[(npts-1)^2 1]));

% Here the physical mode is defined as the mode with the smallest wave
% error, i.e. the mode that most closely matches with the continuous
% eigenvector exp(i \theta \tilde{x})
ix_phys = zeros(length(thetax)*(npts-1),1);

% Isolate the physical mode
for j = 1:length(thetax)
    jevals = evals(:,j);
    jevecs = evecs(:,:,j);

    for k = 1:npts-1
        if (norm(jevecs(:,k)) == 0)
            error('Invalid eigenvector');
        end
    end

    for t = 1:(npts-1)
    for s = 1:(npts-1)

        % Exact mode with wavenumber (s,t)
        % (in 2D there are (npts-1)^2 wavenumbers)
        sx = (s-1)*(npts-1) + t;
        thetax_s = thetax(j) + 2 * pi * (s-1);
        if (thetax_s > pi*(npts-1))
            thetax_s = thetax_s - 2*pi*(npts-1);
        end
        thetay_t = thetay(j) + 2 * pi * (t-1);
        if (thetay_t > pi*(npts-1))
            thetay_t = thetay_t - 2*pi*(npts-1);
        end
        exactev = reshape(exp(1i * thetax_s * gll(1:npts-1)) * exp(1i * thetay_t * gll(1:npts-1).'), [1 (npts-1)^2]);
        exactev = exactev / sqrt(abs(exactev*W*exactev'));

        % There are three wave modes associated with this wavenumber; one
        % geostrophically balanced mode (m=1), one right-propagating mode
        % (m=2), and one left-propagating mode (m=3).
        for m = 1:3
            mx = sx + (m-1)*(npts-1)^2;

            % Find the mode with maximum agreement with this wave mode
            snormmax = 0;
            nmax = 0;
            for n = 1:3*(npts-1)^2

                if (norm(jevecs(:,n)) == 0)
                    continue;
                end
                if ((m == 1) && (abs(imag(jevals(n))) > 1.0e-10))
                    continue;
                end
                if ((m == 2) && (imag(jevals(n)) <= -1.0e-10))
                    continue;
                end
                if ((m == 3) && (imag(jevals(n)) >= 1.0e-10))
                    continue;
                end

                % Calculate the level of agreement between the exact
                % eigenvector and the eigenvector produced by the analysis
                nonzeromodes = 0;
                snorm = 0;
                jevecs_h = jevecs(1:(npts-1)^2,n);
                jevecs_u = jevecs((npts-1)^2+1:2*(npts-1)^2,n);
                jevecs_v = jevecs(2*(npts-1)^2+1:3*(npts-1)^2,n);

                if (norm(jevecs_h) > 1.0e-10)
                    nonzeromodes = nonzeromodes + 1;
                    jevecs_h = jevecs_h / sqrt(abs(jevecs_h'*W*jevecs_h));
                    snorm = snorm + abs(exactev * W * conj(jevecs_h));
                end
                if (norm(jevecs_u) > 1.0e-10)
                    nonzeromodes = nonzeromodes + 1;
                    jevecs_u = jevecs_u / sqrt(abs(jevecs_u'*W*jevecs_u));
                    snorm = snorm + abs(exactev * W * conj(jevecs_u));
                end
                if (norm(jevecs_v) > 1.0e-10)
                    nonzeromodes = nonzeromodes + 1;
                    jevecs_v = jevecs_v / sqrt(abs(jevecs_v'*W*jevecs_v));
                    snorm = snorm + abs(exactev * W * conj(jevecs_v));
                end
                if (nonzeromodes == 0)
                    error('logic error: zero eigenvector found');
                end

                snorm = snorm / nonzeromodes;
                if (snorm > snormmax)
                    snormmax = snorm;
                    nmax = n;
                end
            end

            % Store mode of maximum agreement and remove from contention
            if (nmax ~= 0)
                jevecs(:,nmax) = 0;
                pevals(m,t,(s-1)*length(thetax)+j) = jevals(nmax);
            end
        end
    end
    end
end
pevals;
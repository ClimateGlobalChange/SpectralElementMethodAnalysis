function [evals, evecs] = se_isolate_physical_mode(npts, theta, evals, evecs)
% se_isolate_physical_mode - Sort the eigenvalues in order of increasing
% physical wavenumber, effectively mapping the wave modes to their
% physically relevant counterparts.
%
% Syntax:  [evals, evecs] = se_isolate_physical_mode(npts, theta, evals, evecs)
%
% Inputs:
%    npts - Order of the spectral element method
%    theta - Vector of non-dimensional wavenumbers (k \Delta x_e)
%    evals - Set of eigenvalues (npts-1 x length(theta))
%    evecs - Set of eigenvectors (npts-1 x npts-1 x length(theta))
%
% Remarks:
%    This function is normally called from se_eig.
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

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[gll,gwt] = lglnodes(npts,0,1);
W = diag(gwt(1:npts-1));

% Here the physical mode is defined as the mode with the smallest wave
% error, i.e. the mode that most closely matches with the continuous
% eigenvector exp(i \theta \tilde{x})
ix_phys = zeros(length(theta)*(npts-1),1);

for j = 1:length(theta)
    jevals = evals(:,j);
    jevecs = evecs(:,:,j);

    for k = 1:npts-1
        if (norm(jevecs(:,k)) == 0)
            error('Invalid eigenvector');
        end
    end

    mx = 0;
    for m = 1:npts-1
        if (mx == 0)
            dir = 1;
            thetax = (theta(j) + 2 * pi * (m-1));
            if (thetax > (npts-1)*pi)
                mx = m;
            end
        end
        if (mx ~= 0)
            dir = -1;
            thetax = (theta(j) - 2 * pi * (m-mx+1));
        end
        exactev = exp(1i * thetax * gll(1:npts-1));

        if (dir == -1)
            dir;
        end
        allwaveerr = zeros(npts-1,1);
        for k = 1:npts-1
            if (norm(jevecs(:,k)) == 0)
                allwaveerr(k) = 20.0;
                continue;
            end
            discrev = jevecs(:,k).';
            allwaveerr(k) = 1 - abs(discrev*W*conj(exactev)) / sqrt(abs(discrev*W*discrev')) / sqrt(abs(exactev'*W*exactev));
            if (imag(jevals(k)) * dir < -1e-10)
                allwaveerr(k) = 10.0;
            end
        end

        for k = 1:npts-1
            if (allwaveerr(k) == min(allwaveerr))
                if (dir == 1)
                    evals(m,j) = jevals(k);
                    evecs(:,m,j) = jevecs(:,k);
                else
                    evals(npts-1-(m-mx),j) = jevals(k);
                    evecs(:,npts-1-(m-mx),j) = jevecs(:,k);
                end
                jevecs(:,k) = 0;
                break;
            end
        end
    end
end
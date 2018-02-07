function T = time_discrete_M(time_scheme, M, dt)
% time_discrete_M - Calculate the discrete evolution matrix for the
% specified temporal integration scheme.  That is, for semi-discrete
% equation dphi_j/dt = M phi_j, compute matrix T for discrete update
% (phi_j)^(n+1) = T (phi_j)^n.
%
% Syntax:  T = time_discrete_M(time_scheme, M, dt)
%
% Inputs:
%    time_scheme - Temporal integration rule to apply (see below)
%    M - Semi-discrete update matrix (np x np)
%    dt - Discrete timestep size
%
% Outputs:
%    T - Discrete evolution matrix (np x np)
%
% time_scheme values allowed:
%    'FE'      Forward Euler
%    'RK2'     2nd-order Runge-Kutta
%    'RK3'     3rd-order Runge-Kutta
%    'RK4'     4th-order Runge-Kutta
%    'SSPRK32' Strong stability preserving 3-stage 2nd-order Runge-Kutta
%    'SSPRK43' Strong stability preserving 4-stage 3rd-order Runge-Kutta
%    'SSPRK53' Strong stability preserving 5-stage 3rd-order Runge-Kutta
%    'SSPRK63' Strong stability preserving 6-stage 3rd-order Runge-Kutta
%    'SSPRK73' Strong stability preserving 7-stage 3rd-order Runge-Kutta
%    'SSPRK54' Strong stability preserving 5-stage 4th-order Runge-Kutta
%    'KG53'    5-stage 3rd-order Kinnmark and Gray
%
% Example:
%    >> theta = pi/2;
%    >> [M1,M2,M3] = se_derivative_matrices(3);
%    >> M = M1 * exp(-1i * theta) + M2 + M3 * exp(1i * theta);
%    >> T = time_discrete_M('RK3', M, 0.5)
%
%    T =
%       0.3750 - 0.3125i   1.0417 + 0.5417i
%      -0.2708 + 0.5208i   0.5000 + 0.0833i
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 08-Oct-2017

%------------- BEGIN CODE --------------

if (length(size(M)) ~= 2)
    error('Invalid semi-discrete update matrix M');
end
if (size(M,1) ~= size(M,2))
    error('Invalid semi-discrete update matrix M');
end

npts = size(M,1);

if (strcmp(time_scheme,'FE'))
    T = eye(npts) + dt * M; 
elseif (strcmp(time_scheme,'RK2'))
    T = eye(npts) + dt * M + dt^2/2 * M^2;
elseif (strcmp(time_scheme,'RK3'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3;
elseif (strcmp(time_scheme,'SSPRK32'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3;
elseif (strcmp(time_scheme,'RK4'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 + dt^4/24 * M^4;
elseif (strcmp(time_scheme,'RK5'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 ...
        + dt^4/24 * M^4 + dt^5/120 * M^5;
elseif (strcmp(time_scheme,'SSPRK43'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 ...
        + dt^4 / 48 * M^4;
elseif (strcmp(time_scheme,'SSPRK53'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 ...
        + 0.03143907625 * dt^4 * M^4 ...
        - 0.002372197238 * dt^5 * M^5;
elseif (strcmp(time_scheme,'SSPRK63'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 ...
        + 0.03512028968 * dt^4 * M^4 ...
        - 0.003992765626 * dt^5 * M^5 ...
        + 0.0001891377878 * dt^6 * M^6;
elseif (strcmp(time_scheme,'SSPRK73'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 ...
        + 0.03745121297 * dt^4 * M^4 ... 
        - 0.005240485244 * dt^5 * M^5 ...
        + 0.0004073846037 * dt^6 * M^6 ...
        - 0.00001357253392 * dt^7 * M^7;
elseif (strcmp(time_scheme,'SSPRK54'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 + dt^4/24 * M^4 ...
        - 0.004477718301 * dt^5 * M^5;
elseif (strcmp(time_scheme,'KG53'))
    T = eye(npts) + dt * M + dt^2/2 * M^2 + dt^3/6 * M^3 ...
        + dt^4/30 * M^4 ...
        + dt^5/150 * M^5;
    
else
    error('Invalid temporal scheme');
end

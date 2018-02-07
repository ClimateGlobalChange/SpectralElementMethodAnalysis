function [MDx,MDy] = se2d_deriv_matrices(npts)
% se2d_deriv_matrices - Generate matrices for the spectral element
% first derivative operator in the x and y coordinate directions.
%
% Syntax:  [MDx,MDy] = se2d_deriv_matrices(npts)
%
% Inputs:
%    npts - Order of the spectral element method
%
% Outputs:
%    MDx - Matrix of coefficients for the x derivative operator for all
%          neighboring spectral elements (npts-1 x npts-1 x 3 x 3)
%    MDy - Matrix of coefficients for the y derivative operator for all
%          neighboring spectral elements (npts-1 x npts-1 x 3 x 3)
%
% Remarks:
%    Degrees of freedom within the 2D spectral element are first given in
%    the x direction, then in the y direction. 
%
% Example usage:
%    >> [MDx,MDy] = se2d_deriv_matrices(4);
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 02-Feb-2018

%------------- BEGIN CODE --------------
if (~isscalar(npts))
    error('npts argument must be a scalar');
end

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[myroots,myweights] = lglnodes(npts,0,1);

% Find coefficients of characteristic polynomials
[legpoly,difflegpoly] = characteristic_polyfit(myroots);

% Obtain mass matrix
M = zeros(npts*npts,1);
Mast = zeros(npts*npts,1);
for n = 1:npts
for m = 1:npts
    nbasis = (m-1)*npts+n;
    M(nbasis) = myweights(n) * myweights(m);
    Mast(nbasis) = M(nbasis);

    if ((n == 1) || (n == npts))
        M(nbasis) = M(nbasis) * 2;
    end
    if ((m == 1) || (m == npts))
        M(nbasis) = M(nbasis) * 2;
    end
end
end

% Obtain differentiation matrix
Dx = zeros(npts*npts,npts*npts);
Dy = zeros(npts*npts,npts*npts);
for n = 1:npts
for m = 1:npts
    nbasis = (m-1)*npts+n;

    for s = 1:npts
    for t = 1:npts
        mbasis = (t-1)*npts+s;

        if (m == t)
            for p = 1:npts
                Dx(nbasis,mbasis) = Dx(nbasis,mbasis) ...
                    + myweights(p) * myweights(m) ...
                        * polyval(difflegpoly(n,:), myroots(p)) ...
                        * polyval(legpoly(s,:), myroots(p));
            end
        end
        if (n == s)
            for q = 1:npts
                Dy(nbasis,mbasis) = Dy(nbasis,mbasis) ...
                    + myweights(n) * myweights(q) ...
                        * polyval(difflegpoly(m,:), myroots(q)) ...
                        * polyval(legpoly(t,:), myroots(q));
            end
        end
    end
    end
end
end

MDx = se_extract_submatrices(M, Dx, npts);
MDy = se_extract_submatrices(M, Dy, npts);

myweights;

end

% Extract the submatrices 
function MDall = se_extract_submatrices(M, D, npts)

Dast = D;

D(1:npts,1:npts) = D(1:npts,1:npts) + Dast((npts-1)*npts+1:npts*npts,(npts-1)*npts+1:npts*npts);
D((npts-1)*npts+1:npts*npts,(npts-1)*npts+1:npts*npts) = D((npts-1)*npts+1:npts*npts,(npts-1)*npts+1:npts*npts) + Dast(1:npts,1:npts);

% Mirror integrals along edges of constant x
D(1:npts:npts*npts,1:npts:npts*npts) = D(1:npts:npts*npts,1:npts:npts*npts) + Dast(npts:npts:npts*npts,npts:npts:npts*npts);
D(npts:npts:npts*npts,npts:npts:npts*npts) = D(npts:npts:npts*npts,npts:npts:npts*npts) + Dast(1:npts:npts*npts,1:npts:npts*npts);

% Mirror integrals diagonally
D(1,1) = D(1,1) + Dast(npts*npts,npts*npts);
D(npts*npts,npts*npts) = D(npts*npts,npts*npts) + Dast(1,1);
D(npts,npts) = D(npts,npts) + Dast((npts-1)*npts+1,(npts-1)*npts+1);
D((npts-1)*npts+1,(npts-1)*npts+1) = D((npts-1)*npts+1,(npts-1)*npts+1) + Dast(npts,npts);

% Obtain product
MD = inv(diag(M)) * D;


% Extract unidirectional (in Y) submatrices
MD1 = MD(1:(npts-1)*npts,1:(npts-1)*npts);

MD0 = zeros(size(MD1));
MD0(1:npts,:) = MD((npts-1)*npts+1:npts*npts,1:(npts-1)*npts);

MD2 = zeros(size(MD1));
MD2(:,1:npts) = MD(1:(npts-1)*npts,(npts-1)*npts+1:npts*npts);

% Extract all submatrices
MD00 = zeros((npts-1)*(npts-1));
MD01 = zeros((npts-1)*(npts-1));
MD02 = zeros((npts-1)*(npts-1));
MD10 = zeros((npts-1)*(npts-1));
MD11 = zeros((npts-1)*(npts-1));
MD12 = zeros((npts-1)*(npts-1));
MD20 = zeros((npts-1)*(npts-1));
MD21 = zeros((npts-1)*(npts-1));
MD22 = zeros((npts-1)*(npts-1));

ix1 = [];
for n = 1:npts-1
    ix1 = [ix1 (n-1)*npts+1:n*npts-1];
end
ix0 = [];
for n = 1:npts-1
    ix0 = [ix0 n*npts];
end

MD00(ix0-npts+1,:) = MD0(ix0,ix1);
MD01 = MD0(ix1,ix1);
MD02(:,ix0-npts+1) = MD0(ix1,ix0);

MD10(1:(npts-1):(npts-1)*(npts-1),:) = MD1(ix0,ix1);
MD11 = MD1(ix1,ix1);
MD12(:,1:(npts-1):(npts-1)*(npts-1)) = MD1(ix1,ix0);

MD20(ix0-npts+1,:) = MD2(ix0,ix1);
MD21 = MD2(ix1,ix1);
MD22(:,ix0-npts+1) = MD2(ix1,ix0);

MDall = zeros((npts-1)*(npts-1),(npts-1)*(npts-1),3,3);
MDall(:,:,1,1) = MD00;
MDall(:,:,1,2) = MD01;
MDall(:,:,1,3) = MD02;
MDall(:,:,2,1) = MD10;
MDall(:,:,2,2) = MD11;
MDall(:,:,2,3) = MD12;
MDall(:,:,3,1) = MD20;
MDall(:,:,3,2) = MD21;
MDall(:,:,3,3) = MD22;

end
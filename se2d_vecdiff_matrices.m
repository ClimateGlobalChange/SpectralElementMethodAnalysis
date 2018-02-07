function [DDall_uu, DDall_uv, DDall_vu, DDall_vv, VDall_uu, VDall_uv, VDall_vu, VDall_vv] = se2d_vecdiff_matrices(npts)
% se2d_deriv_matrices - Generate matrices for the spectral element
% divergence damping and vorticity damping operators.
%
% Syntax:  [DDall_uu, DDall_uv, DDall_vu, DDall_vv, VDall_uu, VDall_uv, VDall_vu, VDall_vv] = se2d_vecdiff_matrices(npts)
%
% Inputs:
%    npts - Order of the spectral element method
%
% Outputs:
%    DDall_uu - Matrix of coefficients for the divergence damping operator
%          affecting the u field in terms of the u field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%    DDall_uv - Matrix of coefficients for the divergence damping operator
%          affecting the u field in terms of the v field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%    DDall_vu - Matrix of coefficients for the divergence damping operator
%          affecting the v field in terms of the u field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%    DDall_vv - Matrix of coefficients for the divergence damping operator
%          affecting the v field in terms of the v field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%    VDall_uu - Matrix of coefficients for the vorticity damping operator
%          affecting the u field in terms of the u field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%    VDall_uv - Matrix of coefficients for the vorticity damping operator
%          affecting the u field in terms of the v field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%    VDall_vu - Matrix of coefficients for the vorticity damping operator
%          affecting the v field in terms of the u field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%    VDall_vv - Matrix of coefficients for the vorticity damping operator
%          affecting the v field in terms of the v field for all
%          neigbhoring spectral elements (npts-1 x npts-1 x 3 x 3)
%
% Remarks:
%    Degrees of freedom within the 2D spectral element are first given in
%    the x direction, then in the y direction. 
%
% Example usage:
%    >> [DDall_uu, DDall_uv, DDall_vu, DDall_vv, VDall_uu, VDall_uv,
%    VDall_vu, VDall_vv] = se2d_vecdiff_matrices(4);
%
% Author: Paul Ullrich
% University of California, Davis
% Email address: paullrich@ucdavis.edu
% Last revision: 02-Feb-2018

%------------- BEGIN CODE --------------
if (~isscalar(npts))
    error('npts argument must be a scalar');
end
if (~exist('xonly'))
    xonly = 0;
end
if ((xonly ~= 0) && (xonly ~= 1))
    error('invalid value for agument xonly; value should be 0 or 1');
end

% Compute Gauss-Legendre-Lobatto points on [0,1] interval
[myroots,myweights] = lglnodes(npts,0,1);

% Find coefficients of characteristic polynomials
[~,difflegpoly] = characteristic_polyfit(myroots);

% Obtain mass matrix
M = zeros(npts*npts,1);
for n = 1:npts
for m = 1:npts
    nbasis = (m-1)*npts+n;
    M(nbasis) = myweights(n) * myweights(m);

    if ((n == 1) || (n == npts))
        M(nbasis) = M(nbasis) * 2;
    end
    if ((m == 1) || (m == npts))
        M(nbasis) = M(nbasis) * 2;
    end
end
end

% Construct open curl and divergence operator
Curl = zeros(npts*npts,2*npts*npts);
Div = zeros(npts*npts,2*npts*npts);
for n = 1:npts
for m = 1:npts
    nbasis = (m-1)*npts+n;

    for s = 1:npts
        mbasis = (m-1)*npts+s;

        % DxU
        Div(nbasis,mbasis) = polyval(difflegpoly(s,:), myroots(n));

        % - DxU
        Curl(nbasis,npts*npts+mbasis) = polyval(difflegpoly(s,:), myroots(n));

        mbasis = (s-1)*npts+n;
        
        % DyV
        Div(nbasis,npts*npts+mbasis) = polyval(difflegpoly(s,:), myroots(m));

        % DyU
        Curl(nbasis,mbasis) = - polyval(difflegpoly(s,:), myroots(m));
    end
end
end

% Construct closed gradient operator
Grad = zeros(2*npts*npts,npts*npts);
GradT = zeros(2*npts*npts,npts*npts);

for n = 1:npts
for m = 1:npts
    nbasis = (m-1)*npts+n;

    for s = 1:npts
        % Dx component of grad
        mbasis = (m-1)*npts+s;
        Grad(nbasis,mbasis) = - polyval(difflegpoly(n,:), myroots(s)) * myweights(s) / myweights(n);

        % Dy component of grad
        mbasis = (s-1)*npts+n;
        Grad(npts*npts+nbasis,mbasis) = - polyval(difflegpoly(m,:), myroots(s)) * myweights(s) / myweights(m);
    end
    
end
end

GradT(1:npts*npts,:) = - Grad(npts*npts+1:2*npts*npts,:);
GradT(npts*npts+1:2*npts*npts,:) = Grad(1:npts*npts,:);

% Construct divergence damping operator
DDlocal = Grad * Div;

% Construct vorticity damping operator
VDlocal = GradT * Curl;

% Apply DSS
DD = apply_DSS_vector(npts, DDlocal);
VD = apply_DSS_vector(npts, VDlocal);

% Extract submatrices
npx = npts*npts;

DDall_uu = extract_submatrices(npts, DD(1:npx,1:npx));
DDall_uv = extract_submatrices(npts, DD(1:npx,npx+1:2*npx));
DDall_vu = extract_submatrices(npts, DD(npx+1:2*npx,1:npx));
DDall_vv = extract_submatrices(npts, DD(npx+1:2*npx,npx+1:2*npx));

VDall_uu = extract_submatrices(npts, VD(1:npx,1:npx));
VDall_uv = extract_submatrices(npts, VD(1:npx,npx+1:2*npx));
VDall_vu = extract_submatrices(npts, VD(npx+1:2*npx,1:npx));
VDall_vv = extract_submatrices(npts, VD(npx+1:2*npx,npx+1:2*npx));

myweights;
end

% Apply a virtual DSS operation to the vector operator DD
function DD = apply_DSS_vector(npts, DD)
npx = npts*npts;
DD(1:npx,1:npx) = apply_DSS(npts, DD(1:npx,1:npx));
DD(1:npx,npx+1:2*npx) = apply_DSS(npts, DD(1:npx,npx+1:2*npx));
DD(npx+1:2*npx,1:npx) = apply_DSS(npts, DD(npx+1:2*npx,1:npx));
DD(npx+1:2*npx,npx+1:2*npx) = apply_DSS(npts, DD(npx+1:2*npx,npx+1:2*npx));
end

% Apply a virtual DSS operation to the operator D
function D = apply_DSS(npts, D)

Dast = D;

% Mirror integrals along edges of constant y
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

% Reweight
D(1:npts,1:npts*npts) = D(1:npts,1:npts*npts) / 2;
D((npts-1)*npts+1:npts*npts,1:npts*npts) = D((npts-1)*npts+1:npts*npts,1:npts*npts) / 2;
D(1:npts:npts*npts,1:npts*npts) = D(1:npts:npts*npts,1:npts*npts) / 2;
D(npts:npts:npts*npts,1:npts*npts) = D(npts:npts:npts*npts,1:npts*npts) / 2;

end

% Extract the submatrices of DD needed for linear analysis
function MDall = extract_submatrices(npts, MD)

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
MD02(:,1:(npts-1):(npts-1)*(npts-1)) = MD0(ix1,ix0);

MD10(1:(npts-1):(npts-1)*(npts-1),:) = MD1(ix0,ix1);
MD11 = MD1(ix1,ix1);
MD12(:,1:(npts-1):(npts-1)*(npts-1)) = MD1(ix1,ix0);

MD20(1:(npts-1):(npts-1)*(npts-1),:) = MD2(ix0,ix1);
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
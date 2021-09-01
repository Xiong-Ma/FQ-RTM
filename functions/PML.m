function [pmlz, pmlzh, pmlx, pmlxh] = PML(npml, free, omega, nz, nx, dh)
% PML: calculate the reciprocal of PML coefficients
%
% [pmlz, pmlzh, pmlx, pmlxh] = PML(npml, free, omega, nz, nx, dh)
%
% Input :
%       npml  = the thickness of PML [gridpoints]
%       free  = the top boundry is free-surface boundary(=1), or pml boundary(=0)
%       omega = the angular frequency
%       nz    = the number of grid-point in the z-direction
%       nx    = the number of grid-point in the x-direction
%       dh    = the spatial sampling interval [m], dh=dx=dz
% Output :
%       pmlz  = the reciprocal of PML coefficients in the z-direction
%       pmlzh = the reciprocal of half PML coefficients in the z-direction
%       pmlx  = the reciprocal of PML coefficients in the x-direction
%       pmlxh = the reciprocal of half PML coefficients in the x-direction
%

a0_pml=1.79; % the constant in PML coefficients calculation

% define thickness of PML [m]
lPML = npml*dh;

% calculate PML damping profiles sigma_x in the x-direction
sigma_x = zeros(nx,1);
sx      = zeros(nx,1);
for i=1:nx
    sigma_x(i) = 0.0;
    % define damping profile at left PML boundary
    if(i <= npml)
        sigma_x(i) = omega.*a0_pml.*((npml-i).*dh./lPML).^2;
%         sigma_x(i) = omega.*a0_pml.*((npml-i+1).*dh./lPML).^2;
    end
    % define damping profile at right PML boundary
    if(i >= nx - npml)
        sigma_x(i) = omega.*a0_pml.*((dh.*(nx-npml+1-i))./lPML).^2;
%         sigma_x(i) = omega.*a0_pml.*((dh.*(nx-npml-i))./lPML).^2;
    end
    sx(i) = 1-(sigma_x(i)*1i/omega);
end

% calculate PML damping profiles sigma_z in the z-direction
sigma_z = zeros(nz,1);
sz      = zeros(nz,1);
for j=1:nz
    sigma_z(j) = 0.0;
    % define damping profile at top PML boundary
    if((free==0)&&(j <= npml))
        sigma_z(j) = omega.*a0_pml.*((npml-j).*dh./lPML).^2;
%         sigma_z(j) = omega.*a0_pml.*((npml-j+1).*dh./lPML).^2;
    end
    % define damping profile at bottom PML boundary
    if(j >= nz-npml)
        sigma_z(j) = omega.*a0_pml.*((dh.*(nz-npml+1-j))./lPML).^2;
%         sigma_z(j) = omega.*a0_pml.*((dh.*(nz-npml-j))./lPML).^2;
    end
    sz(j) = 1-(sigma_z(j)*1i/omega);
end

% the PML coefficients in the half-point grid
zh = sz;
for j=2:nz-1
    zh(j) = (sz(j) + sz(j+1))./2.0;
end
xh = sx;
for i=2:nx-1
    xh(i) = (sx(i) + sx(i+1))./2.0;
end

% the reciprocal of PML coefficients
pmlz=1./sz;
pmlx=1./sx;
pmlzh=1./zh;
pmlxh=1./xh;

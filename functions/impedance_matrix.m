function A = impedance_matrix(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v, Q, fr)
% IMPEDANCE_MATRIX: generate a impedance matrix use an optimal 9-point, finite-differnce method
%
% A = impedance_matrix(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v)
% A = impedance_matrix(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v, Q)
% A = impedance_matrix(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v, Q, fr)
%
% Input :
%       flag  = 1 for acoustic, 2 for dissipation, 3 for dispersion, 4 for viscoacoustic
%       nz    = the number of grid-point in the z-direction
%       nx    = the number of grid-point in the x-direction
%       dh    = the spatial sampling interval [m], dh=dx=dz
%       pmlz  = the reciprocal of PML coefficients in the z-direction
%       pmlzh = the reciprocal of half PML coefficients in the z-direction
%       pmlx  = the reciprocal of PML coefficients in the x-direction
%       pmlxh = the reciprocal of half PML coefficients in the x-direction
%       fk    = the current frequency silce to be calculated [Hz]
%       v     = the velocity model
%       Q     = the Q model
%       fr    = the reference frequency used for modeling the viscoacoustic effects [Hz]
%               [default fr = 50 Hz]
% Output :
%       A = the spare impedance matrix (stroes in a compressed format)
%
% Please refer:
% Jo, C. H., C. Shin, and J. H. Suh, 1996, An optimal 9-point, finite-differnce, frequency-space,
%        2-D scalar wave extrapolator: Geophysics, 61(2), 529¨C537.
% Bernhard Hustedt, Stephane Operto and Jean Virieux, 2004, Mixed-grid and staggered-grid finite-difference
%        methods for frequency-domain acoustic wave modelling: Geophys. J. Int, 157, 1269¨C1296
%
% The 9-point we used are Northwest(NW) point, West(W) point, Southwest(SW) point,
%                         North(N) point, Central(C) point, South(S) point,
%                         Northeast(NE) point, East(E) point, Southeast(SE) point,
%                             NW       N      NE
%                               * ---- * ---- *
%                               |    C |      |
%                             W * ---- * ---- * E
%                               |      |      |
%                               * ---- * ---- *
%                            SW        S      SE
%
% by Xiong Ma and Bo Gao, April 2021
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
%
%
if(nargin<12)
   fr = 50;
end
if(nargin<11)
   Q = v; 
end
if (fk==0) % 0 Hz component is not calculate in the frequency domain method
    disp('0 Hz component is not calculate, exit !');
    return;
end

% the square of the inverse complex velocity
iv2 = (1./v).^2;
switch flag
    case 1
        iv2 = (1./v).^2;   % acoustic
    case 2
        iv2 = ((1./v).*(1-1i./2./Q)).^2;   % dissipation-only
    case 3
        iv2 = ((1./v).*(1-(1./pi./Q).*log(fk/fr))).^2;  % dispersion-only
    case 4
        iv2 = ((1./v).*(1-(1./pi./Q).*log(fk/fr)).*(1-1i./2./Q)).^2; % viscoacoustic
end

% the optimal coefficient
a = 0.5461;
c = 0.6248;
d = 0.09381;
e = (1-c-4*d)/4;

% auxiliary variables
idh2   = 1.0/dh^2;
coff1  = 0.5*(1+a)*idh2;
coff2  = 0.25*(1-a)*idh2;
omega2 = (2.0*pi*fk)^2;
kq2    = omega2.*iv2;

% Assemble impedance matrix
% loop over Cartesian grid
k = 1;
count = 1;
row_ind = zeros(9*nx*nz,1);
col_ind = zeros(9*nx*nz,1);
value_A = zeros(9*nx*nz,1);

for i=1:nx
    for j=1:nz
        % implement discrete FDFD operator ...
        % NW gridpoint
        if((i > 1) && (j > 1))
            lump = pmlz(j)*pmlzh(j-1) + pmlx(i)*pmlxh(i-1);
            value_A(count) = e*kq2(j-1,i-1) + coff2*lump;
            row_ind(count) = k;
            col_ind(count) = k-1-nz;
            count = count + 1;
        end
        
        % W gridpoint
        if(i > 1)
            lump1 = pmlx(i)*pmlxh(i-1);
            if(j == 1)
                lump2 = pmlz(j)*pmlzh(j);
            else
                lump2 = pmlz(j)*pmlzh(j-1) + pmlz(j)*pmlzh(j);
            end
            value_A(count) = d*kq2(j,i-1) + coff1*lump1 - coff2*lump2;
            row_ind(count) = k;
            col_ind(count) = k-nz;
            count = count + 1;
        end
        
        % SW gridpoint
        if((i > 1) && (j < nz))
            lump = pmlz(j)*pmlzh(j) + pmlx(i)*pmlxh(i-1);
            value_A(count) = e*kq2(j+1,i-1) + coff2*lump;
            row_ind(count) = k;
            col_ind(count) = k+1-nz;
            count = count + 1;
        end
        
        % N gridpoint
        if(j > 1)
            lump1 = pmlz(j)*pmlzh(j-1);
            if(i == 1)
                lump2 = pmlx(i)*pmlxh(i);
            else
                lump2 = pmlx(i)*pmlxh(i-1) + pmlx(i)*pmlxh(i);
            end
            value_A(count) = d*kq2(j-1,i) + coff1*lump1 - coff2*lump2;
            row_ind(count) = k;
            col_ind(count) = k-1;
            count = count + 1;
        end
        
        % Central gridpoint
        if ((i > 1) && (j > 1))
            lump = pmlz(j)*pmlzh(j-1) + pmlz(j)*pmlzh(j) + pmlx(i)*pmlxh(i-1) + pmlx(i)*pmlxh(i);
        end
        if ((i > 1) && (j == 1))
            lump = pmlz(j)*pmlzh(j) + pmlx(i)*pmlxh(i-1) + pmlx(i)*pmlxh(i);
        end
        if ((i == 1) && (j > 1))
            lump = pmlz(j)*pmlzh(j-1) + pmlz(j)*pmlzh(j) + pmlx(i)*pmlxh(i);
        end
        if ((i == 1) && (j == 1))
            lump = pmlz(j)*pmlzh(j) + pmlx(i)*pmlxh(i);
        end
        value_A(count) = c*kq2(j,i) - coff1*lump;
        row_ind(count) = k;
        col_ind(count) = k;
        count = count + 1;
        
        % S gridpoint
        if (j < nz)
            lump1 = pmlz(j)*pmlzh(j);
            if(i == 1)
                lump2 = pmlx(i)*pmlxh(i);
            else
                lump2 = pmlx(i)*pmlxh(i-1) + pmlx(i)*pmlxh(i);
            end
            value_A(count) = d*kq2(j+1,i) + coff1*lump1 - coff2*lump2;
            row_ind(count) = k;
            col_ind(count) = k+1;
            count = count + 1;
        end
        
        % NE gridpoint
        if((j > 1) && (i < nx))
            lump = pmlz(j)*pmlzh(j-1) + pmlx(i)*pmlxh(i);
            value_A(count) = e*kq2(j-1,i+1) + coff2*lump;
            row_ind(count) = k;
            col_ind(count) = k-1+nz;
            count = count + 1;
        end
        
        % E gridpoint
        if(i < nx)
            lump1 = pmlx(i)*pmlxh(i);
            if(j == 1)
                lump2 = pmlz(j)*pmlzh(j);
            else
                lump2 = pmlz(j)*pmlzh(j-1) + pmlz(j)*pmlzh(j);
            end
            value_A(count) = d*kq2(j,i+1) + coff1*lump1 - coff2*lump2;
            row_ind(count) = k;
            col_ind(count) = k+nz;
            count = count + 1;
        end
        
        % SE gridpoint
        if((i < nx) && (j < nz))
            lump = pmlz(j)*pmlzh(j) + pmlx(i)*pmlxh(i);
            value_A(count) = e*kq2(j+1,i+1) + coff2*lump;
            row_ind(count) = k;
            col_ind(count) = k+1+nz;
            count = count + 1;
        end
        k = k + 1;
    end
end

zeros_ind = find(value_A==0);
value_A(zeros_ind) = [];
row_ind(zeros_ind)  = [];
col_ind(zeros_ind)  = [];

% assemble sparse matrix A
A = sparse(row_ind,col_ind,value_A);
end

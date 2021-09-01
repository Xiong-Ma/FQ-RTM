function [nz, nx, v_pml, Q_pml] = model_pml(nz, nx, npml, v, Q)
% MODEL_PML: Extend model parameters with PML boundary
%
% [nz, nx, v_pml, Q_pml] = model_pml(nz, nx, npml, v, Q)
%
% Input :
%       nz   = the number of grid-point in the z-direction
%       nx   = the number of grid-point in the x-direction
%       npml = the thickness of PML [gridpoints]
%       v    = the original velocity model
%       Q    = the original Q model
% Output :
%       nz    = the number of grid-point in the z-direction of the extended model
%       nx    = the number of grid-point in the x-direction of the extended model
%       v_pml = the original velocity model with PML boundary
%       Q_pml = the original Q model with PML boundary
%

% the number of grid-point with PML boundary
nx_new = nx + 2 * npml;
nz_new = nz + 2 * npml;

% the velocity, Q models with PML boundary
v_pml = zeros(nz_new, nx_new);
Q_pml = zeros(nz_new, nx_new);

% extend velocity model
% central model (equal to the original model)
v_pml(1 + npml : nz+npml, 1+npml : nx+npml) = v;
% left boundary (equal to left boundry of the original model)
for i=1:npml
    v_pml(:,i) = v_pml(:,1+npml);
end
% right boundary (equal to right boundry of the original model)
for i=(nx + npml + 1):nx_new
    v_pml(:,i) = v_pml(:,nx + npml);
end
% top boundary (equal to top boundry of the original model)
for j=1:npml
    v_pml(j,:) = v_pml(npml + 1,:);
end
% bottom boundary
for j=(nz + npml + 1):nz_new
    v_pml(j,:) = v_pml(nz + npml,:);
end

% extend Q model
% central model (equal to the original model)
Q_pml(1+npml:nz+npml,1+npml:nx+npml) = Q;
% left boundary (equal to left boundry of the original model)
for i=1:npml
    Q_pml(:,i) = Q_pml(:,1+npml);
end
% right boundary
for i=(nx + npml + 1):nx_new
    Q_pml(:,i) = Q_pml(:,nx + npml);
end
% top boundary
for j=1:npml
    Q_pml(j,:) = Q_pml(npml + 1,:);
end
% bottom boundary
for j=(nz + npml + 1):nz_new
    Q_pml(j,:) = Q_pml(nz + npml,:);
end

nx = nx_new;
nz = nz_new;
end

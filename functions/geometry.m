function [nshots, shot_index, nrecs, rec_index] = geometry(nz, nx, npml)
% Define acquisition geometry
%
% [nshots, shot_index, nrecs, rec_index] = geometry(nz, nx, npml)
%
% Input :
%       nz   = the number of grid-point in the z-direction
%       nx   = the number of grid-point in the x-direction
%       npml = the thickness of PML [gridpoints]
% Output :
%       nshots     = the shot number
%       shot_index = the index of shots
%       nrecs      = the receiver number
%       shot_index = the index of receivers
%

% define the intervel of shots and receivers
  s_int = 10;  % the shots intervel
  r_int = 1;   % the receivers intervel

% calculate discrete shot and receiver positions on FD grid
  shot_nx = 1:s_int:nx;
  shot_nz = 1;
  rec_nx  = 1:r_int:nx;
  rec_nz  = 1;
  
% the number of shots and receivers  
  nshots = length(shot_nx);
  nrecs  = length(rec_nx);
    
% calculate the index of shots and receivers on 1D FD grid
  shot_index = zeros(nshots,1);
  for i = 1:nshots
      shot_index(i) = (shot_nx(i)+npml-1)*(nz+2*npml)+shot_nz+npml;
  end
  rec_index = zeros(nrecs,1);
  for i=1:nrecs
      rec_index(i)=(rec_nx(i)+npml-1)*(nz+2*npml)+rec_nz+npml;
  end
end

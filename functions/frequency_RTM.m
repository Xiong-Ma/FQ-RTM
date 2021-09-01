function image = frequency_RTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources, v)
% FREQUENCY_RTM: the frequency domain reverse-time migration (RTM) of prestack shot gathers (depth migration)
%
% image = frequency_RTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources, v)
%
% The migration can be divided into three parts. Firstly, calculate the forward-propagating of the source wavefield
% Then, compute the reverse-time-propagated of the receiver wavefield. Finally, apply the cross-correction imaging
% condition of two wavefields to obatin the migrated profile.
%
% Input :
%       nz    = the number of grid-point in the z-direction
%       nx    = the number of grid-point in the x-direction
%       dh    = the spatial sampling interval [m], dh=dx=dz
%       npml  = the thickness of PML [gridpoints]
%       pmlz  = the reciprocal of PML coefficients in the z-direction
%       pmlzh = the reciprocal of half PML coefficients in the z-direction
%       pmlx  = the reciprocal of PML coefficients in the x-direction
%       pmlxh = the reciprocal of half PML coefficients in the x-direction
%       w_f   = the source wavelet in the frequency domain, can be a vector or scalar,
%       f     = the frequency vecter or scalar, if w_f is a vector, it should be a vector with the same length
%       shot_index = the shot index in 1D grid
%       rec_index = the receiver index in 1D grid (default: rec_index=shot_index)
%       receiver_sources  = the frequency-domain prestack seismic data as virtual sources at receivers' location,
%                           it has three demension (1--frequency, 2--receiver, 3--shot)
%       v     = the velocity model (with PML)
% Output :
%       image = the image result
%
%

display('Beginning: Frequency domain reverse-time migration......');
tic;

% check the parameters
[nzv, nxv] = size(v);
if(nzv~=nz || nxv~=nx)
    disp('The model size is not equal to the parameters nz and/or nx!');
    return;
end
clear nzv nxv;

npmlz = length(pmlz);
npmlzh = length(pmlzh);
npmlx = length(pmlx);
npmlxh = length(pmlxh);
if(npmlz~=npmlzh || npmlz~=nz || npmlx~=npmlxh || npmlx~=nx)
    disp('The length of PML coefficients is not equal to the parameters nz and/or nx!');
    return;
end
clear npmlz npmlzh npmlx npmlxh;

nw = length(w_f);
nf = length(f);
if(nw~=nf)
    disp('Warnning: The length of wavelet and frequency is not equal!');
    return;
end
clear nw nf;

% the parameters used for allocating memory
num_grid = nz*nx; % the number of grid point
nshots   = length(shot_index); % the number of shots
nres     = length(rec_index);  % the number of receivers

% check the data and parameters
if( length(receiver_sources(1,:,1)) ~= nres)
    disp('The number of receivers (prestack data) is not equal to the geometry (forward modeling parameter)!');
    return;
end
if( length(receiver_sources(1,1,:)) ~= nshots)
    disp('The number of shot gathers (prestack data) is not equal to the geometry (forward modeling parameter)!');
    return;
end

% the 0 Hz component is not calculate, the wavefield=0
nonzero_ind=find(f~=0);

% the cross-correction imaging condition
image = zeros(num_grid,1);   % the image of full source and full receiver wavefields
tmp   = zeros(num_grid,1);   % the source illumination

loopNum = length(nonzero_ind);
parfor_progress(loopNum);  % Initialize

parfor k=nonzero_ind
    fk  = f(k);  % the current frequency
    
    % generate impedance matrix at current frequency
    A = impedance_matrix(1, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v);
    [L,U,P,O]=lu(A);         % LU decomposition for sparse matrix
    
    % calculate the frequency domain forward-propagating source wavefield
    wfk = -w_f(k);    % the source wavelet's spectrum at current frequency
    S = sparse(shot_index,1:nshots,wfk,num_grid,nshots); % the shot term (the right hand term of the equation Ax=S)
    S = full(S);
    source_wavefield=(O*(U\(L\(P*S)))); % the elimination and backsubstitution steps for solving source wavefield

    % calculate the frequency domain backward-propagating receiver wavefield
    S(rec_index,:) = -receiver_sources(k,:,:); % the receiver sources term (the right hand term of the equation Ax=R)
    receiver_wavefield = (O*(U\(L\(P*S)))); % the elimination and backsubstitution steps for solving receiver wavefield
    
    % apply the normalized cross-correction imaging condition 
    image = image + sum((-2*pi*fk)*source_wavefield.*receiver_wavefield, 2);
    tmp = tmp + sum(source_wavefield.*conj(source_wavefield), 2); % source illumination
    
    parfor_progress;  % Count
end

image = image./tmp;    % source illumination (scale the image profile)
image = real(image);
image = reshape(image,nz,nx);
h = fspecial('laplacian',0.0);  % Laplacian filter
image=imfilter(image,h,'replicate'); % Laplacian filtering
image = image(npml+1:nz-npml,npml+1:nx-npml);

parfor_progress(0);  % Clean up  

toc;
display('Finished: Frequency domain reverse-time migration......');
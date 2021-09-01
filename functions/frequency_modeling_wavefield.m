function freq_wavefield_full = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v, Q, fr)
% FREQUENCY_MODELING_WAVEFIELD: the frequency domain forward modeling for generating the mutli-shots seismic wavefields
%
% freq_wavefield_full = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v)
% freq_wavefield_full = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v, Q)
% freq_wavefield_full = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v, Q, fr)
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
%       w_f   = the source wavelet in the frequency domain, can be a vector or scalar, 
%       f     = the frequency vecter or scalar, if w_f is a vector, it should be a vector with the same length         
%       shot_index = the shot index in 1D grid
%       v     = the velocity model
%       Q     = the Q model
%       fr    = the reference frequency used for modeling the viscoacoustic effects [Hz]
%               [default fr = 50 Hz]
% Output :
%       freq_wavefield_full = the frequency-domain wavefield which is a 3D dataset,
%                       1st dimension is frequency, 2nd dimension is wavefield, 3rd dimension is shot indix  
%
% Note that this function is just suitable for small model because we record all frequency wavefields,
% so it requires a lot of memory. Usually, we only record the wavefields at the receivers location.
% Please refer the function "frequency_modeling_record" which only record the wavefields at the receivers location.
%
% by Xiong Ma and Bo Gao, April 2021
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
%
%

display('Beginning: Frequency domain forward modeling for generating seismic wavefields......');
tic;

if(nargin<14)
   fr = 50;
end
if(nargin<13)
   Q = v; 
end

% the 0 Hz component is not calculate, the wavefield=0
nonzero_ind=find(f~=0); 

% the parameters used for allocating memory
num_grid = nz*nx;   % the number of grid point
nfs      = length(f);           % the length of f 
nshots   = length(shot_index);  % the number of shots

freq_wavefield_full=zeros(nfs,num_grid,nshots);  % the frequency domain wavefields

loopNum = length(nonzero_ind);
parfor_progress(loopNum);  % Initialize
parfor k=nonzero_ind
    fk = f(k);  % the current frequency
    wfk= -w_f(k);% the wavelet's spectrum at current frequency
    
    % generate impedance matrix at current frequency
    A = impedance_matrix(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v, Q, fr); 
    [L,U,P,O]=lu(A);         % LU decomposition for sparse matrics
    
    % the sources term (the right hand term of the equation Ax=S)
    S = sparse(shot_index,1:nshots,wfk,num_grid,nshots); 
    S = full(S);
    
    S = O*(U\(L\(P*S))); % the elimination and backsubstitution steps of LU solution (wavefields)
    
    freq_wavefield_full(k,:,:)=S; % record the wavefields (iterate for all frequency components)
    
    parfor_progress;  % Count
end

parfor_progress(0);  % Clean up  

toc;
display('Finished: Frequency domain forward modeling for generating seismic wavefields......');
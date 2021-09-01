function image = frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources, v, Q, fr, sigma)
% FREQUENCY_RTM: the frequency-domain Q-compensation reverse-time migration (RTM) of prestack shot gathers (depth migration)
%
% image = frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources, v, Q, sigma)
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
%       v     = the velocity model
%       Q     = the Q model
%       fr    = the reference frequency used for modeling the viscoacoustic effects [Hz]
%               [default fr = 50 Hz]
%       sigma = the stabilization factor, [default sigma = 10^(-5)]
% Output :
%       image = the image result
%

display('Beginning: Frequency domain Q-compensation reverse-time migration......');
tic;

if(nargin<17)
   sigma = 10^(-5);
end
if(nargin<16)
   fr = 50;
end

% check the parameters
[nzv, nxv] = size(v);
[nzq, nxq] = size(Q);
if(nzv~=nzq || nzv~=nz || nxv~=nxq || nxv~=nx)
    disp('The model size is not equal to the parameters nz and/or nx!');
    return;
end
clear nzv nxv nzq nxq;

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
nshots   = length(shot_index);    % the number of shots
nres     = length(rec_index);     % the number of receivers

if( length(receiver_sources(1,:,1)) ~= nres)
    disp('The number of receivers (prestack data) is not equal to the geometry (forward modeling parameter)!');
    return;
end
if( length(receiver_sources(1,1,:)) ~= nshots)
    disp('The number of shot gathers (prestack data) is not equal to the geometry (forward modeling parameter)!');
    return;
end

% the 0 Hz component is not calculate, the wavefield=0
nonzero_ind = find(f~=0); 

% the cross-correction imaging condition
image = zeros(num_grid,1);   % the image
tmp   = zeros(num_grid,1);   % the source illumination

gauss_filter = fspecial('gaussian',3,1); % the gaussian smooth filter

loopNum = length(nonzero_ind);
parfor_progress(loopNum);  % Initialize

parfor k=nonzero_ind
    fk  = f(k);  % the current frequency
    
    % generate impedance matrix for dispersion-only modeling
    A = impedance_matrix(3, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v, Q, fr);
    [LP,UP,PP,OP]=lu(A); % LU decomposition for sparse matrix (dispersion-only modeling)
    
    % generate impedance matrix for viscoacoustic modeling
    A = impedance_matrix(4, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, fk, v, Q, fr); 
    [LQ,UQ,PQ,OQ]=lu(A); % LU decomposition for sparse matrix (viscoacoustic modeling)
    
    % calculate the frequency domain backward-propagating receiver wavefield
    S = zeros(num_grid, nshots);  % initialize the receiver source term for all shots (the right hand term)
    S(rec_index,:) = -receiver_sources(k,:,:); % the receiver sources at the current frequency
    % the elimination and backsubstitution steps for calculation the backward-propagating receiver wavefield
    wavefield_dispr = (OP*(UP\(LP\(PP*S)))); % dispersion-only source wavefield
    wavefield_visco = (OQ*(UQ\(LQ\(PQ*S)))); % viscoacoustic source wavefield
    % calculate the stabilised source compensation wavefield
    a = abs(wavefield_visco).*abs(wavefield_dispr);
    b = abs(wavefield_visco).*abs(wavefield_visco);
    wavefield_ampli = zeros(num_grid,nshots); % initialize the stabilised amplitude compensation operator
    for i=1:nshots
        aa = reshape(a(:,i), nz, nx);
        bb = reshape(b(:,i), nz, nx);
        a(:,i) = reshape(imfilter(aa,gauss_filter), nz*nx, 1); % use the gaussian smooth filter to smooth the wavefield
        b(:,i) = reshape(imfilter(bb,gauss_filter), nz*nx, 1);
        wavefield_ampli(:,i)  = a(:,i)./(b(:,i)+max(b(:,i))*sigma); % calculate the stabilised amplitude compensation operator
    end
    % the stabilised receiver compensation wavefield (amplitude+phase)
    receiver_wavefield = wavefield_ampli.*wavefield_dispr; 

    % calculate the frequency domain forward-propagating source wavefield
    wfk = -w_f(k);% the wavelet's spectrum at current frequency
    S = sparse(shot_index, 1:nshots, wfk, num_grid, nshots); % the shot term (the right hand term of the equation Ax=S)
    S = full(S);
    % the elimination and backsubstitution steps for calculation the forward-propagating source wavefield
    wavefield_dispr = (OP*(UP\(LP\(PP*S))));   % dispersion-only source wavefield
    wavefield_visco = (OQ*(UQ\(LQ\(PQ*S))));   % viscoacoustic source wavefield
    % calculate the stabilised source compensation wavefield
    a = abs(wavefield_visco).*abs(wavefield_dispr);
    b = abs(wavefield_visco).*abs(wavefield_visco);
    for i=1:nshots
        aa = reshape(a(:,i), nz, nx);
        bb = reshape(b(:,i), nz, nx);
        a(:,i) = reshape(imfilter(aa,gauss_filter), nz*nx, 1); % use the gaussian smooth filter to smooth the wavefield
        b(:,i) = reshape(imfilter(bb,gauss_filter), nz*nx, 1);
        wavefield_ampli(:,i)  = a(:,i)./(b(:,i)+max(b(:,i))*sigma);  % the stabilised amplitude compensation wavefield
    end
    % the stabilised source compensation wavefield (amplitude+phase)
    source_wavefield = wavefield_ampli.*wavefield_dispr; 

    % apply the normalized cross-correction imaging condition
    image = image + sum((-2*pi*fk)*source_wavefield.*receiver_wavefield, 2);
    tmp   = tmp + sum(wavefield_dispr.*conj(wavefield_dispr), 2); % use for scaling the image profile

    parfor_progress;  % Count
end
image=image./tmp; % scale the image profile
image = real(image);
image = reshape(image, nz, nx);
h = fspecial('laplacian',0.0);   % Laplacian filter
image=imfilter(image,h,'replicate'); % Laplacian filtering
image = image(npml+1:nz-npml, npml+1:nx-npml);

parfor_progress(0);  % Clean up  

toc;
display('Finished: Frequency domain Q-compensation reverse-time migration......');
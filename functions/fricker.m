function [w_f, f]=fricker(dt, wf, wl)
% FRICKER: creates a causal Ricker wavelet in the frequncy domain 
%
% [w_f, f]=fricker(dt, wf, wl)
% [w_f, f]=fricker(dt, wf)
% [w_f, f]=fricker(dt) 
%
% Input :
%       dt = the temporal sample rate [s] 
%       wf = the wavelet's peak frequency [Hz] (default: 30 Hz)
%       wl = the wavelet's length [sample-point] (default: 501)
% Output :
%       w_f = the frequncy domain Ricker wavelet
%       f = the frequency vecter
%
% Notes:
% The amplitude of the Ricker wavelet at a frequency of 2.5*fpeak is 
% approximately 4 percent of that at the dominant frequency fpeak.
% The Ricker wavelet effectively begins at time t = -1.0/fpeak.  Therefore,
% for practical purposes, a causal wavelet may be obtained by a time delay
% of 1.0/fpeak.
% The Ricker wavelet has the shape of the second derivative of a Gaussian.
%

if(nargin<3)
   wl=501;
 end
 if(nargin<2)
   wf=30; 
 end
 
%  time delay
  td=1/wf;
  
 % create a time vector
  tmax=dt*(wl-1);
  tw= 0.:dt:tmax;
  tw=tw-td;

 % create the time domain wavelet
  x = pi*wf.*tw;
  xx = x.*x;
  wavelet=(1-2*xx).*exp(-xx);
  
% create a frequency vector
  fnyq=1./(2*(tw(2)-tw(1)));   % the Nyquist frequency 
  f=linspace(0.,fnyq,length(tw)/2+1);
  
% the frequency domain wavelet
  fft_w = fft(wavelet,wl);
  w_f=fft_w(1:length(f));
  
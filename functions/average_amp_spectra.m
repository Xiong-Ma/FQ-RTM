function [f, amp] = average_amp_spectra(d, dt, nf, fremax)
% avg_amp_spectra: the average amplitude spectrum of a 2D plofile
%
% [f, amp] = average_amp_spectra(flag, d, dt, nf, fremax)
%
% Input :
%       d     = 2D seismic data
%       dt    = the time-stepping interval [s]
%       nf    = the length of the Fourier transform
%       fremax= the maximum frequency displayed in the spectrum
%
% Output :
%       f    = the frequency vector
%       amp  = the average amplitude spectrum of the 2D seismic data
%

[~,N] = size(d);
D = fft(d,nf,1);
df = 1/(nf*dt);
nfremax = floor(fremax/df)+1;
f = (0:nfremax-1)*df;
f = f';
amp = sum(abs(D(1:nfremax,:)),2)/N;

end
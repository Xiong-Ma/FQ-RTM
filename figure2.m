clear all;clc;close all;
addpath('./functions','./parfor_progress');
%% The model parameters
dt = 0.002;   % the time sampling interval [s]
nt = 501;     % the time record points
fpeak = 30;      % the modeling (peak) frequency [Hz]
omega = 2.0*pi*fpeak;  % the angular frequency
fr = 200;     % the reference frequency for modeling the viscoacoustic effects

% define the modeling grid
nz = 200;   % the number of grid-point in the z-direction
nx = 200;   % the number of grid-point in the x-direction
dh = 3;    % the spatial sampling interval [m], dh=dx=dz

% define the parameter modeles (velocity, Q)
v = ones(nz,nx)*1000;  % the velocity model
Q = ones(nz,nx)*10;    % the Q model

% the parameters of PML absorption boundary
npml = 20;  % the thickness of PML [gridpoints]
free = 0;  % the top boundry is free-surface boundary(=1), or pml boundary(=0)

% the shots and receivers locations
shot_nx = 101;  % the shot positions on the modeling grid
shot_nz = 101;
shot_index=zeros(length(shot_nx),1); % the shot index in the 1D model vector
for i=1:length(shot_nx)
    shot_index(i)=(shot_nx(i)+npml-1)*(nz+2*npml)+shot_nz+npml;
end
  
% define the source wavelet
[w_f, f]=fricker(dt, fpeak, nt); % a causal Ricker wavelet in the frequncy domain

% extend the original model to contain the PML boundary
[nz, nx, v, Q] = model_pml(nz, nx, npml, v, Q); % add pml bounary to the original model 
[pmlz, pmlzh, pmlx, pmlxh] = PML(npml, free, omega, nz, nx, dh); % calculate the absorption coefficients

%% 
% To improve efficiency, we select limited frequency band for modeling rather than the full frequency band.
% We select a high-cut frequency and the frequency higher than it will be ignored
fhigh = 101;
for i = 1:length(f)
    if(fhigh-f(i)<dt)
        break;
    end
end
f = f(1:i);
w_f = w_f(1:i);

%% The forward modeling for acoustic equation in the frequency domain
display('############# The forward modeling for acoustic equation #############');
flag = 1; 
freq_wavefields_acoustic = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v);

%% The forward modeling for dissipation-only equation in the frequency domain
display('############# The forward modeling for dissipation-only equation #############');
flag = 2;  
freq_wavefields_dissipation = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v, Q, fr);

%% The forward modeling for for dispersion-only equation in the frequency domain
display('############# The forward modeling for dispersion-only equation #############');
flag = 3;  
freq_wavefields_dispersion = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v, Q, fr);

%% The forward modeling for viscoacoustic equation in the frequency domain
display('############# The forward modeling for viscoacoustic equation #############');
flag = 4;  
freq_wavefields_viscoacoustic = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v, Q, fr);

%% transform the frequency domain wavefields to time domain (snapshots)
time_wavefields_acoustic      = ifft(freq_wavefields_acoustic,nt,1,'symmetric');
time_wavefields_dissipation   = ifft(freq_wavefields_dissipation,nt,1,'symmetric');
time_wavefields_dispersion    = ifft(freq_wavefields_dispersion,nt,1,'symmetric');
time_wavefields_viscoacoustic = ifft(freq_wavefields_viscoacoustic,nt,1,'symmetric');

%% Choose t=280ms (dt=0.002, that is, tk=141) for comparsion
tk=141;

time_wavefield_acoustic_tk = time_wavefields_acoustic(tk,:,1);
time_wavefield_acoustic_tk_2D = reshape(time_wavefield_acoustic_tk,nz,nx);
clear time_wavefield_acoustic_tk;

time_wavefield_dissipation_tk = time_wavefields_dissipation(tk,:,1);
time_wavefield_dissipation_tk_2D=reshape(time_wavefield_dissipation_tk,nz,nx);
clear time_wavefield_dissipation_tk;

time_wavefield_dispersion_tk = time_wavefields_dispersion(tk,:,1);
time_wavefield_dispersion_tk_2D = reshape(time_wavefield_dispersion_tk,nz,nx);
clear time_wavefield_dispersion_tk;

time_wavefield_viscoacoustic_tk = time_wavefields_viscoacoustic(tk,:,1);
time_wavefield_viscoacoustic_tk_2D = reshape(time_wavefield_viscoacoustic_tk,nz,nx);
clear time_wavefield_viscoacoustic_tk;

% compare the snapshots in four cases 
a1 = time_wavefield_acoustic_tk_2D(npml+1:(nz/2),npml+1:(nx/2));
a2 = time_wavefield_dissipation_tk_2D(npml+1:(nz/2),npml+1:(nx/2));
a2 = imrotate(a2,-90);
a3 = time_wavefield_dispersion_tk_2D(npml+1:(nz/2),npml+1:(nx/2));
a3 = imrotate(a3,90);
a4 = time_wavefield_viscoacoustic_tk_2D(npml+1:(nz/2),npml+1:(nx/2));
a4 = imrotate(a4,180);
aa=[a1,a2;a3,a4];
[nz,nx]=size(aa);
z=1:dh:dh*nz;
x=1:dh:dh*nx;

alpha=0:pi/360:2*pi;
xx=300+265*cos(alpha);
yy=300+265*sin(alpha);

yyy=ones(size(x))*300;
xxx=ones(size(z))*300;

figure;
imagesc(x,z,aa);
colormap(gray);
hold on;
plot(xx,yy,'r--');
hold on;
scatter(300,300,'r*')
hold on;
plot(x,yyy,'k-.',xxx,z,'k-.');
set(gca,'xtick',0:100:dh*nx,'fontsize',12);
set(gca,'ytick',0:100:dh*nz,'fontsize',12); 
ylabel('Distance (m)');
xlabel('Distance (m)');
set(gcf,'Position',[100 100 350 320]); 
% set(gca,'position',[0.75 0.1 0.2 0.81]);
caxis([-.1,.1]);
strings1={'a) Acoustic'};
annotation('textbox',[0.17,0.88,0.3,0.05],'LineStyle','none',...
   'LineWidth',1,'String',strings1,'fontsize',11);
strings1={'b) Dissipation'};
annotation('textbox',[0.62,0.88,0.3,0.05],'LineStyle','none',...
   'LineWidth',1,'String',strings1,'fontsize',11);
strings1={'c) Dispersion'};
annotation('textbox',[0.17,0.165,0.3,0.05],'LineStyle','none',...
   'LineWidth',1,'String',strings1,'fontsize',11);
strings1={'d) Viscoacoustic'};
annotation('textbox',[0.57,0.165,0.35,0.05],'LineStyle','none',...
   'LineWidth',1,'String',strings1,'fontsize',11);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.5 3.2]);
print ./figures/figure2.eps -depsc2 -r600;
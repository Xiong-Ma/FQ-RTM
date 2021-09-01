%% Important Things: 
%  This script is just to make the reader more aware of the process of frequency-domain forward modeling.
%  In fact, the frequency-domain forward modeling is more efficient for multi-shot survey.
%  In addition, the frequency-domain full wavefields are the intermediate variables.
%  Actually, we do not need to store full wavefields in memory or disk, 
%  but in this test script we store them for the convenient to show wavefields at each frequency.
%  
clear all;clc;close all;
addpath('./functions','./parfor_progress');
%% load and show the velocity model and Q model
load ./data/layered_v.mat;  % the velocity model
load ./data/layered_Q.mat;   % the Q model
[nz,nx]=size(v);

% define the modeling grid
dh = 5;    % the spatial sampling interval [m], dh=dx=dz
z = dh.*(0:nz-1)/1000;
x = dh.*(0:nx-1)/1000;

figure;
imagesc(x,z,v/1000);
colormap(gray);
colorbar('location','Eastoutside');
colorbar;
ylabel(colorbar,'Velocity (km/s)','Fontsize',13);
set(gcf,'Position',[100 100 600 400]); 
set(gca,'Position',[.1 .03 .72 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',13); 
set(gca,'ytick',0:0.2:2,'fontsize',13); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',13);
ylabel('Depth (km)', 'fontsize',13);
set(gca,'xaxislocation','top');

% store the figure in the folder named 'Fig'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 4.0]);
print ./figures/figure5a.eps -depsc2 -r600;

figure;
imagesc(x,z,Q);
colorbar;
colormap(gray);
ylabel(colorbar,'Qulity factor','Fontsize',13);
set(gcf,'Position',[100 100 600 400]); 
set(gca,'Position',[.1 .03 .72 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',13); 
set(gca,'ytick',0:0.2:2,'fontsize',13); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',13);
ylabel('Depth (km)', 'fontsize',13);
set(gca,'xaxislocation','top');

% store the figure in the folder named 'Fig'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 4.0]);
print ./figures/figure5b.eps -depsc2 -r600;

%% The model parameters
dt = 0.002;   % the time sampling interval [s]
nt = 501;    % the time record length [sample points]
t = 0:dt:(nt-1)*dt;
fpeak = 30;       % the modeling (peak) frequency [Hz]
omega = 2.0*pi*fpeak;  % the angular frequency
fr = 100;     % the reference frequency for modeling the viscoacoustic effects

% the parameters of PML absorption boundary
npml = 20;  % the thickness of PML [gridpoints]
free = 0;

% Define the source wavelet
[w_f, f]=fricker(dt, fpeak, nt); % a causal Ricker wavelet in the frequncy domain

%% Define the geometry (the shots and receivers locations)
shot_nx = 101;  % the shot positions on FD grid
shot_nz = 1;
nshots = 1;
shot_index = (shot_nx+npml-1)*(nz+2*npml) + shot_nz + npml; % the shot index in the 1D model vector

r_int = 2;   % the receivers intervel
rec_nx  = 1:r_int:nx;  % the receivers positions on FD grid
rec_nz  = 1;
nrecs  = length(rec_nx);
rec_index = zeros(nrecs,1);
for i=1:nrecs
    rec_index(i)=(rec_nx(i)+npml-1)*(nz+2*npml)+rec_nz+npml;
end

%% extend the original model to contain the PML boundary
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

%% frequency-domain forward modeling for acoustic equation
display('############# The forward modeling for acoustic equation #############');
flag = 1;  % for acoustic equation
freq_wavefields_acoustic = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v);

%% frequency-domain forward modeling for viscoacoustic equation
display('############# The forward modeling for viscoacoustic equation #############');
flag = 4;  % for viscoacoustic equation
freq_wavefields_viscoacoustic = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, v, Q, fr);

%% make a comparison of two wavefields at 50 Hz
fk = 51;
freq_wavefield_acou_fk = freq_wavefields_acoustic(fk,:);
freq_wavefield_acou_fk_2D = reshape(freq_wavefield_acou_fk,nz,nx);
clear freq_wavefield_acou_fk;

freq_wavefield_viscoaco_fk = freq_wavefields_viscoacoustic(fk,:);
freq_wavefield_viscoaco_fk_2D = reshape(freq_wavefield_viscoaco_fk,nz,nx);
clear freq_wavefield_viscoaco_fk;

% compare the real part of two wavefields
figure;
subplot(1,2,1);
imagesc(x,z,real(freq_wavefield_acou_fk_2D(npml+1:nz-npml,npml+1:nx-npml)));
set(gcf,'Position',[100 100 700 300]); 
set(gca,'Position',[.08 .03 .4 .82]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-5,5]);
strings={'a)'};
annotation('textbox',[.01,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);
subplot(1,2,2);
imagesc(x,z,real(freq_wavefield_viscoaco_fk_2D(npml+1:nz-npml,npml+1:nx-npml)));
set(gca,'Position',[.57 .03 .4 .82]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-5,5]);
strings={'b)'};
annotation('textbox',[.5,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7.0 3.0]);
print ./figures/figure6.eps -depsc2 -r600;

% compare the amplitude component of two wavefields
figure;
subplot(1,2,1);
imagesc(x,z,abs(freq_wavefield_acou_fk_2D(npml+1:nz-npml,npml+1:nx-npml)));
set(gcf,'Position',[100 100 700 300]); 
set(gca,'Position',[.08 .03 .4 .82]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([1,5]);
strings={'a)'};
annotation('textbox',[.01,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);
subplot(1,2,2);
imagesc(x,z,abs(freq_wavefield_viscoaco_fk_2D(npml+1:nz-npml,npml+1:nx-npml)));
set(gca,'Position',[.57 .03 .4 .82]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([1,5]);
strings={'b)'};
annotation('textbox',[.5,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7.0 3.0]);
print ./figures/figure7.eps -depsc2 -r600;

%% transform the frequency-domain wavefields to time-domain domain
time_wavefields_acoustic = ifft(freq_wavefields_acoustic,nt,1,'symmetric');% acoustic wavefield
time_wavefields_viscoacoustic = ifft(freq_wavefields_viscoacoustic,nt,1,'symmetric');% viscoacoustic wavefield

%% show the snapshot at t=280ms
tk = 141;
time_wavefield_acou_tk = reshape(time_wavefields_acoustic(tk,:),nz,nx);
time_wavefield_viscoaco_tk = reshape(time_wavefields_viscoacoustic(tk,:),nz,nx);

figure;
subplot(1,2,1);
imagesc(x,z,time_wavefield_acou_tk(npml+1:nz-npml,npml+1:nx-npml));
colormap(gray);
set(gcf,'Position',[100 100 700 300]); 
set(gca,'Position',[.08 .03 .4 .82]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-0.5,1]);
strings={'a)'};
annotation('textbox',[.01,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);
subplot(1,2,2);
imagesc(x,z,time_wavefield_viscoaco_tk(npml+1:nz-npml,npml+1:nx-npml));
colormap(gray);
set(gca,'Position',[.57 .03 .4 .82]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-0.5,1]);
strings={'b)'};
annotation('textbox',[.5,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7.0 3.0]);
print ./figures/figure8.eps -depsc2 -r600;

%% cut-off the direct wave using auxiliary velocity model
va=v(1,1).*ones(size(v));  % auxiliary velocity model

display('############# Cut-off the direct wave for acoustic equation #############');
% frequency domain acoustic wavefield after eliminating the direct wave
flag = 1;  % for acoustic equation
freq_wavefields_acoustic_aux = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, va);
freq_wavefields_acoustic = freq_wavefields_acoustic - freq_wavefields_acoustic_aux;
clear freq_wavefields_acoustic_aux;

display('############# Cut-off the direct wave for viscoacoustic equation #############');
% frequency domain viscoacoustic wavefield after eliminating the direct wave
flag = 4;  % for viscoacoustic equation
freq_wavefield_viscoacoustic_aux = frequency_modeling_wavefield(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, va, Q, fr);
freq_wavefields_viscoacoustic = freq_wavefields_viscoacoustic - freq_wavefield_viscoacoustic_aux;
clear freq_wavefield_viscoacoustic_aux;

%% frequency-domain acoustic record and viscoacoustic record
freq_record_acoustic = freq_wavefields_acoustic(:,rec_index);
freq_record_viscoaco = freq_wavefields_viscoacoustic(:,rec_index);

%% transform the frequency domain records to time domain
time_record_acoustic = ifft(freq_record_acoustic,nt,1,'symmetric');  % acoustic shot gathers
time_record_viscoaco = ifft(freq_record_viscoaco,nt,1,'symmetric');  % viscoacoustic shot gathers

%% show the shot gathers at X=500 m
figure;
subplot(1,2,1);
imagesc(x,t,time_record_acoustic);
colormap(gray);
set(gcf,'Position',[100 100 550 400]); 
set(gca,'Position',[.09 .03 .39 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Time (s)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-.15 .2]);
strings={'a)'};
annotation('textbox',[.01,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);
subplot(1,2,2);
imagesc(x,t,time_record_viscoaco);
colormap(gray);
set(gca,'Position',[.59 .03 .39 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',0:0.2:2,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Time (s)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-.15 .2]);
strings={'b)'};
annotation('textbox',[.5,0.96,0.05,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',12);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 4.0]);
print ./figures/figure9.eps -depsc2 -r600;

%% Plot the traces located at X=200 m, X=500 m, and X=800 m
X = 21;
trace_acoustic_x200 = time_record_acoustic(:,X);
trace_viscoaco_x200 = time_record_viscoaco(:,X);
X = 51;
trace_acoustic_x500 = time_record_acoustic(:,X);
trace_viscoaco_x500 = time_record_viscoaco(:,X);
X = 81;
trace_acoustic_x800 = time_record_acoustic(:,X);
trace_viscoaco_x800 = time_record_viscoaco(:,X);

figure;
subplot(3,1,1);
plot(t,trace_acoustic_x200,'r','LineWidth',1.5);
hold on;
plot(t,trace_viscoaco_x200,'b--','LineWidth',1.5);
set(gcf,'Position',[100 100 600 500]); 
% set(gca,'Position',[.1 .03 .79 .45]);
set(gca,'xtick',0:0.2:2,'fontsize',11); 
set(gca,'ytick',-0.2:0.1:0.2,'fontsize',11); 
ylabel('Amplitude', 'fontsize',11);
xlabel('Time (s)', 'fontsize',11);
legend('Acoustic','Viscoacoustic','Orientation','Horizontal','Location',[0.62 0.935 0.1 0.05]);
strings={'a)'};
annotation('textbox',[0.02,0.94,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',11);
subplot(3,1,2);
plot(t,trace_acoustic_x500,'r','LineWidth',1.5);
hold on;
plot(t,trace_viscoaco_x500,'b--','LineWidth',1.5);
set(gca,'xtick',0:0.2:2,'fontsize',11); 
set(gca,'ytick',-0.2:0.1:0.2,'fontsize',11); 
ylabel('Amplitude', 'fontsize',11);
xlabel('Time (s)', 'fontsize',11);
strings={'b)'};
annotation('textbox',[0.02,0.64,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',11);
subplot(3,1,3);
plot(t,trace_acoustic_x800,'r','LineWidth',1.5);
hold on;
plot(t,trace_viscoaco_x800,'b--','LineWidth',1.5);
set(gca,'xtick',0:0.2:2,'fontsize',11); 
set(gca,'ytick',-0.2:0.1:0.2,'fontsize',11); 
ylabel('Amplitude', 'fontsize',11);
xlabel('Time (s)', 'fontsize',11);
strings={'c)'};
annotation('textbox',[0.02,0.34,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',11);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 5.0]);
print ./figures/figure10.eps -depsc2 -r600;

clear all;clc;close all;
%%
addpath('./functions','./parfor_progress');
load ./data/gas_chimney_vel.mat; % the velocity model
load ./data/gas_chimney_Q.mat;   % the Q model
[nz,nx]=size(v);
    
% define the modeling aera
dh = 10;    % the spatial sampling interval [m], dh=dx=dz
z = dh.*(0:nz-1)/1000;
x = dh.*(0:nx-1)/1000;

figure;
imagesc(x,z(41:end),v(41:end,:)/1000);
colorbar('location','Eastoutside');
colorbar;
ylabel(colorbar,'Velocity (km/s)','Fontsize',12);
set(gcf,'Position',[100 100 600 300]); 
set(gca,'Position',[.1 .03 .72 .83]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.3);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 3.0]);
print ./figures/figure13a.eps -depsc2 -r600;

figure;
imagesc(x,z(41:end),Q(41:end,:));
colorbar;
ylabel(colorbar,'Qulity factor','Fontsize',12);
set(gcf,'Position',[100 100 600 300]); 
set(gca,'Position',[.1 .03 .72 .83]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.3);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 3.0]);
print ./figures/figure13b.eps -depsc2 -r600;

%% The model parameters
% define the modeling frequency
dt = 0.001;   % the time sampling interval [s]
nt = 2501;    % the time record length [sample points]
t = 0:dt:(nt-1)*dt;
fpeak = 30;       % the modeling (peak) frequency [Hz]
omega = 2.0*pi*fpeak;  % the angular frequency
fr = 100;     % the reference frequency

% the parameters of PML absorption boundary
npml = 20;  % the thickness of PML [gridpoints]
free = 0;

% Define the source wavelet
[w_f, f]=fricker(dt, fpeak, nt);  % a causal Ricker wavelet in the frequncy domain

% Define the geometry
[nshots, shot_index, nrecs, rec_index] = geometry(nz, nx, npml);

% extend the original model to contain the PML boundary
[nz, nx, v, Q] = model_pml(nz, nx, npml, v, Q); % add pml bounary to the original model 
[pmlz, pmlzh, pmlx, pmlxh] = PML(npml, free, omega, nz, nx, dh); % calculate the absorption coefficients

%% To improve efficiency, we select the frequency band for forward modeling.
% % Usually, we select a high-cut frequency, and the frequency higher than it will be ignored
fhigh = 101;
for i = 1:length(f)
    if(fhigh-f(i)<dt)
        break;
    end
end
f = f(1:i);
w_f = w_f(1:i);

%% Frequency-domain forward modeling for acoustic equation
display('############# The forward modeling for acoustic equation #############');
flag = 1;  
freq_record_acoustic = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, v);

% transform the frequency domain records to the time domain
time_record_acoustic = ifft(freq_record_acoustic,nt,1,'symmetric');  % acoustic shot gathers
clear freq_record_acoustic;   % delete the frequency-domain acoustic records

%% Frequency-domain forward modeling for viscoacoustic equation
display('############# The forward modeling for viscoacoustic equation #############');
flag = 4;
freq_record_viscoaco = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, v, Q, fr);

% transform the frequency domain records to the time domain
time_record_viscoaco = ifft(freq_record_viscoaco,nt,1,'symmetric');  % viscoacoustic shot gathers
clear freq_record_viscoaco;   % delete the frequency-domain viscoacoustic records

%% Cut off the direct wave for acoustic equation using the auxiliary velocity model
va=v(1,1).*ones(size(v));  % the auxiliary velocity model

display('############# Cut-off the direct wave for acoustic equation #############');
flag = 1;  % for acoustic equation
freq_record_acoustic_a = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, va);

% transform the frequency domain records to the time domain
time_record_acoustic_a = ifft(freq_record_acoustic_a,nt,1,'symmetric');
clear freq_record_acoustic_a;   % delete the frequency-domain auxiliary records

% cut off the direct wave for acoustic shot gathers
time_record_acoustic = time_record_acoustic - time_record_acoustic_a;
clear time_record_acoustic_a;    % delete the auxiliary records

%% Cut off the direct wave for viscoacoustic equation using the auxiliary velocity model
display('############# Cut-off the direct wave for viscoacoustic equation #############');
flag = 4;  % for viscoacoustic equation
freq_record_viscoaco_a = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, va, Q, fr);

% transform the frequency domain records to the time domain
time_record_viscoaco_a = ifft(freq_record_viscoaco_a,nt,1,'symmetric');
clear freq_record_viscoaco_a;   % delete the frequency-domain auxiliary records

% cut off the direct wave for viscoacoustic shot gathers
time_record_viscoaco = time_record_viscoaco - time_record_viscoaco_a;
clear time_record_viscoaco_a;   % delete the auxiliary records

%% Adding the random noise to the viscoacoustic shot gathers
SNR = 2;
time_record_viscoaco_noise = time_record_viscoaco;
for j = 1:nshots
    time_record_viscoaco_noise(:,:,j) = awgn(time_record_viscoaco(:,:,j),SNR,'measured');
end

%% Show the shot gathers at X=2100 m
shot_indix = 21;
shot_acoustic = time_record_acoustic(:,:,shot_indix);
shot_viscoaco = time_record_viscoaco(:,:,shot_indix);
shot_viscoa_n = time_record_viscoaco_noise(:,:,shot_indix);

figure;
imagesc(x,t,shot_acoustic);
colormap(gray);
set(gcf,'Position',[100 100 300 400]); 
set(gca,'Position',[.17 .03 .77 .85]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Time (s)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-.04 .04]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 4.0]);
print ./figures/figure14a.eps -depsc2 -r600;

figure;
imagesc(x,t,shot_viscoaco);
colormap(gray);
set(gcf,'Position',[100 100 300 400]); 
set(gca,'Position',[.17 .03 .77 .85]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Time (s)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-.04 .04]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 4.0]);
print ./figures/figure14b.eps -depsc2 -r600;

figure;
imagesc(x,t,shot_viscoa_n);
colormap(gray);
set(gcf,'Position',[100 100 300 400]); 
set(gca,'Position',[.17 .03 .77 .85]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Time (s)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-.04 .04]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 4.0]);
print ./figures/figure14c.eps -depsc2 -r600;

%% Show seismic traces at X=2100 m
Xs = 210;
shot_acoustic_xs = shot_acoustic(:,Xs);
shot_viscoaco_xs = shot_viscoaco(:,Xs);
shot_viscoa_n_xs = shot_viscoa_n(:,Xs);

[fs, amp_acoustic_xs] = average_amp_spectra(shot_acoustic, dt, nt, 100);
[~, amp_viscoaco_xs] = average_amp_spectra(shot_viscoaco, dt, nt, 100);
[~, amp_viscoa_n_xs] = average_amp_spectra(shot_viscoa_n, dt, nt, 100);

figure;
subplot(2,1,1);
plot(t(50:end),shot_acoustic_xs(50:end),'r','LineWidth',1.0);
hold on;
plot(t(50:end),shot_viscoaco_xs(50:end),'b--','LineWidth',1.0);
hold on;
plot(t(50:end),shot_viscoa_n_xs(50:end),'g-.','LineWidth',1.0);
set(gcf,'Position',[100 100 600 500]); 
set(gca,'Position',[.12 .6 .85 .34]);
set(gca,'xtick',0:0.4:5,'fontsize',12); 
set(gca,'ytick',-.1:0.05:.1,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Time (s)', 'fontsize',12);
xlim([0.5 2.0]);
ylim([-0.08 0.12]);
legend('Acoustic','Viscoacoustic','Viscoacoustic 2dB',...
      'Orientation','Horizontal');
strings={'a)'};
annotation('textbox',[0.02,0.94,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);
subplot(2,1,2);
plot(fs, 0.5*amp_acoustic_xs,'r','LineWidth',1.5);
hold on;
plot(fs, amp_viscoaco_xs,'b--','LineWidth',1.5);
hold on;
plot(fs, amp_viscoa_n_xs,'g-.','LineWidth',1.5);
set(gca,'Position',[.12 .1 .85 .4]);
xlim([0 100]);
ylim([0 1.7]);
set(gca,'xtick',0:20:200,'fontsize',12); 
set(gca,'ytick',0:0.4:2,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Frequency (Hz)', 'fontsize',12);
legend('Acoustic \times 0.5','Viscoacoustic','Viscoacoustic 2dB',...
      'Orientation','Vertical');
strings={'b)'};
annotation('textbox',[0.02,0.5,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 5.0]);
print ./figures/figure15.eps -depsc2 -r600;

%% Preporcessing for RTM and Q-RTM 
[w_f, f] = fricker(dt, fpeak, nt); % frequncy-domain Ricker wavelet as the forward source 
receiver_sources_a = fft(time_record_acoustic,nt);  % the received acoustic prestack as virtual source
receiver_sources_v = fft(time_record_viscoaco,nt);  % the received viscoacoustic prestack as virtual source
receiver_sources_vn = fft(time_record_viscoaco_noise,nt); % the noisy viscoacoustic prestack as virtual source

receiver_sources_a = conj(receiver_sources_a);
receiver_sources_v = conj(receiver_sources_v);
receiver_sources_vn = conj(receiver_sources_vn);

% % To improve efficiency, we select the frequency band to migration.
% % Usually, we select a high-cut frequency, and the frequency higher than it will be ignored
fhigh = 81;
for i = 1:length(f)
    if(fhigh-f(i)<dt)
        break;
    end
end
f = f(1:1:i);
w_f = w_f(1:1:i);
receiver_sources_a = receiver_sources_a(1:1:i,:,:);
receiver_sources_v = receiver_sources_v(1:1:i,:,:);
receiver_sources_vn = receiver_sources_vn(1:1:i,:,:);

%% Acoustic RTM with acoustic prestack data (reference image)
display('############# FA-RTM with acoustic data #############');
image_ARTM_a=frequency_RTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_a, v);

figure;
imagesc(x,z(41:end),image_ARTM_a(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16a.eps -depsc2 -r600;

%% Acoustic RTM with the noise-free viscoacoustic prestack data
display('############# FA-RTM with viscoacoustic data #############');
image_ARTM_v=frequency_RTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_v, v);

figure;
imagesc(x,z(41:end),image_ARTM_v(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16b.eps -depsc2 -r600;

%% Stabilized Q-compensated RTM with the noise-free viscoacoustic prestack data
display('############# Stabilized FQ-RTM with viscoacoustic data #############');
sigma = 10^(-9);
image_QRTM_v1=frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_v, v, Q, fr, sigma);

figure;
imagesc(x,z(41:end),image_QRTM_v1(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16c.eps -depsc2 -r600;

%% Stabilized Q-compensated RTM with the noise-free viscoacoustic prestack data
display('############# Stabilized FQ-RTM with viscoacoustic data #############');
sigma = 10^(-6);
image_QRTM_v2=frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_v, v, Q, fr, sigma);

figure;
imagesc(x,z(41:end),image_QRTM_v2(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16d.eps -depsc2 -r600;

%% Stabilized Q-compensated RTM with the noise-free viscoacoustic prestack data
display('############# Stabilized FQ-RTM with viscoacoustic data #############');
sigma = 10^(-4);
image_QRTM_v3=frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_v, v, Q, fr, sigma);

figure;
imagesc(x,z(41:end),image_QRTM_v3(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16e.eps -depsc2 -r600;

%% The differences between reference image and Q-compensated images
image_error1 = image_QRTM_v1  - image_ARTM_a;
image_error2 = image_QRTM_v2  - image_ARTM_a;
image_error3 = image_QRTM_v3  - image_ARTM_a;

figure;
imagesc(x,z(41:end),image_error1(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16f.eps -depsc2 -r600;

figure;
imagesc(x,z(41:end),image_error2(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16g.eps -depsc2 -r600;

figure;
imagesc(x,z(41:end),image_error3(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure16h.eps -depsc2 -r600;

%% Show the seismic traces extracted from the images
X = 110;  % x=1100 m
trace_acoustic_x1100 = image_ARTM_a(:,X);
trace_viscoaco_x1100 = image_ARTM_v(:,X);
trace_q_comp_x1100 = image_QRTM_v1(:,X);
X = 210;  % x=2100 m
trace_acoustic_x2100 = image_ARTM_a(:,X);
trace_viscoaco_x2100 = image_ARTM_v(:,X);
trace_q_comp_x2100 = image_QRTM_v1(:,X);
X = 350;  % x=3500 m
trace_acoustic_x3500 = image_ARTM_a(:,X);
trace_viscoaco_x3500 = image_ARTM_v(:,X);
trace_q_comp_x3500 = image_QRTM_v1(:,X);

figure;
subplot(3,1,1);
plot(z(50:end),trace_acoustic_x1100(50:end),'r','LineWidth',1.5);
hold on;
plot(z(50:end),trace_viscoaco_x1100(50:end),'b--','LineWidth',1.5);
hold on;
plot(z(50:end),trace_q_comp_x1100(50:end),'g-.','LineWidth',1.5);
set(gcf,'Position',[100 100 600 500]); 
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',-10000:500:1000,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Depth (km)', 'fontsize',12);
xlim([0.5 2]);
ylim([-500 500]);
legend('Acoustic','Viscoacoustic','Compensated',...
      'Orientation','Horizontal','Location',[0.44 0.93 0.2 0.05]);
strings={'a)'};
annotation('textbox',[0.02,0.94,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);
subplot(3,1,2);
plot(z(50:end),trace_acoustic_x2100(50:end),'r','LineWidth',1.5);
hold on;
plot(z(50:end),trace_viscoaco_x2100(50:end),'b--','LineWidth',1.5);
hold on;
plot(z(50:end),trace_q_comp_x2100(50:end),'g-.','LineWidth',1.5);
xlim([0.5 2]);
ylim([-1300 1300]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',-1000:1000:1000,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Depth (km)', 'fontsize',12);
strings={'b)'};
annotation('textbox',[0.02,0.64,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);
subplot(3,1,3);
plot(z(50:end),trace_acoustic_x3500(50:end),'r','LineWidth',1.5);
hold on;
plot(z(50:end),trace_viscoaco_x3500(50:end),'b--','LineWidth',1.5);
hold on;
plot(z(50:end),trace_q_comp_x3500(50:end),'g-.','LineWidth',1.5);
xlim([0.5 2]);
ylim([-500 300]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',-500:500:500,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Depth (km)', 'fontsize',12);
strings={'c)'};
annotation('textbox',[0.01,0.34,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 5.0]);
print ./figures/figure17.eps -depsc2 -r600;

%% Stabilized Q-compensated RTM with the noisy viscoacoustic prestack data
display('############# Stabilized FQ-RTM with viscoacoustic data #############');
sigma = 10^(-9);
image_QRTM_vn1=frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_vn, v, Q, fr, sigma);

figure;
imagesc(x,z(41:end),image_QRTM_vn1(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure18a.eps -depsc2 -r600;

%% Stabilized Q-compensated RTM with the noisy viscoacoustic prestack data
display('############# Stabilized FQ-RTM with viscoacoustic data #############');
sigma = 10^(-6);
image_QRTM_vn2=frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_vn, v, Q, fr, sigma);

figure;
imagesc(x,z(41:end),image_QRTM_vn2(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure18b.eps -depsc2 -r600;

%% Stabilized Q-compensated RTM with the noisy viscoacoustic prestack data
display('############# Stabilized FQ-RTM with viscoacoustic data #############');
sigma = 10^(-4);
image_QRTM_vn3=frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_vn, v, Q, fr, sigma);

figure;
imagesc(x,z(41:end),image_QRTM_vn3(41:end,:));
colormap(gray);
set(gcf,'Position',[100 100 550 300]); 
set(gca,'Position',[.1 .03 .86 .82]);
set(gca,'xtick',0:0.5:5,'fontsize',12); 
set(gca,'ytick',0:0.5:5,'fontsize',12); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',12);
ylabel('Depth (km)', 'fontsize',12);
set(gca,'xaxislocation','top');
caxis([-500 500]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 3.0]);
print ./figures/figure18c.eps -depsc2 -r600;
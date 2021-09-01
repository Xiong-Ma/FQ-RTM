clear all;clc;close all;
%% load the parameter modeles (velocity, Q)
addpath('./functions','./parfor_progress');
load ./data/layered_v.mat;  % the velocity model
load ./data/layered_Q.mat;   % the Q model
[nz,nx]=size(v);

% define the modeling aera
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

%% define the modeling parameters
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

%% Define the geometry
[nshots, shot_index, nrecs, rec_index] = geometry(nz, nx, npml);

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

%% The frequency-domain forward modeling for acoustic equation
display('############# The forward modeling for acoustic equation #############');
flag = 1;  % for acoustic equation
freq_record_acoustic = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, v);

%% The frequency-domain forward modeling for viscoacoustic equation
display('############# The forward modeling for viscoacoustic equation #############');
flag = 4;  % for viscoacoustic equation
freq_record_viscoaco = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, v, Q, fr);

%% Cut-off the direct wave using auxiliary velocity model
va=v(1,1).*ones(size(v)); % Auxiliary velocity

display('############# Cut-off the direct wave for acoustic equation #############');
flag = 1;  % for acoustic equation
freq_record_acoustic_aux = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, va);
freq_record_acoustic = freq_record_acoustic - freq_record_acoustic_aux; % cut off the direct wave
clear freq_record_acoustic_a;    % delete the auxiliary records

display('############# Cut-off the direct wave for viscoacoustic equation #############');
flag = 4;  % for viscoacoustic equation
freq_record_viscoaco_aux = frequency_modeling_record(flag, nz, nx, dh, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, va, Q, fr);
freq_record_viscoaco = freq_record_viscoaco - freq_record_viscoaco_aux; % cut off the direct wave
clear freq_record_viscoaco_a;   % delete the auxiliary records

%% transform the frequency-domain records to time-domain records
time_record_acoustic = ifft(freq_record_acoustic,nt,1,'symmetric');  % acoustic shot gathers
time_record_viscoaco = ifft(freq_record_viscoaco,nt,1,'symmetric');  % viscoacoustic shot gathers

%% Preporcessing for RTM and Q-RTM in the frequency domain
[w_f, f] = fricker(dt, fpeak, nt); % frequncy-domain Ricker wavelet as the forward source 
receiver_sources_a = fft(time_record_acoustic,nt);  % the received acoustic prestack as virtual source
receiver_sources_v = fft(time_record_viscoaco,nt);  % the received viscoacoustic prestack as virtual source

receiver_sources_a = conj(receiver_sources_a);
receiver_sources_v = conj(receiver_sources_v);

% To improve efficiency, we select the frequency band to migration.
% Usually, we select a high-cut frequency, and the frequency higher than it will be ignored
fhigh = 101;
for i = 1:length(f)
    if(fhigh-f(i)<dt)
        break;
    end
end
f = f(1:i);
w_f = w_f(1:i);
receiver_sources_a = receiver_sources_a(1:i,:,:);
receiver_sources_v = receiver_sources_v(1:i,:,:);

%% Acoustic RTM with acoustic prestack data
display('############# FA-RTM with acoustic data #############');
image_ARTM_acoustic = frequency_RTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_a, v);

figure;
imagesc(x,z,image_ARTM_acoustic);
colormap(gray);
set(gcf,'Position',[100 100 530 400]); 
set(gca,'Position',[.11 .03 .86 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',13); 
set(gca,'ytick',0:0.2:2,'fontsize',13); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',13);
ylabel('Depth (km)', 'fontsize',13);
set(gca,'xaxislocation','top');
caxis([-200 200]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.3 4.0]);
print ./figures/figure11a.eps -depsc2 -r600;

%% Acoustic RTM with viscoacoustic prestack data
display('############# FA-RTM with viscoacoustic data #############');
image_ARTM_viscoac = frequency_RTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_v, v);

figure;
imagesc(x,z,image_ARTM_viscoac);
colormap(gray);
set(gcf,'Position',[100 100 530 400]); 
set(gca,'Position',[.11 .03 .86 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',13); 
set(gca,'ytick',0:0.2:2,'fontsize',13); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',13);
ylabel('Depth (km)', 'fontsize',13);
set(gca,'xaxislocation','top');
caxis([-200 200]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.3 4.0]);
print ./figures/figure11b.eps -depsc2 -r600;

%% Stabilized Q-compensated RTM with viscoacoustic prestack data
display('############# Stabilized FQ-RTM with viscoacoustic data #############');
sigma = 10^(-6);
image_QRTM_viscoac = frequency_QRTM(nz, nx, dh, npml, pmlz, pmlzh, pmlx, pmlxh, w_f, f, shot_index, rec_index, receiver_sources_v, v, Q, fr, sigma);

figure;
imagesc(x,z,image_QRTM_viscoac);
colormap(gray);
set(gcf,'Position',[100 100 530 400]); 
set(gca,'Position',[.11 .03 .86 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',13); 
set(gca,'ytick',0:0.2:2,'fontsize',13); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',13);
ylabel('Depth (km)', 'fontsize',13);
set(gca,'xaxislocation','top');
caxis([-200 200]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.3 4.0]);
print ./figures/figure11c.eps -depsc2 -r600;

%% The differences between reference image and Q-compensated image
image_error = image_QRTM_viscoac - image_ARTM_acoustic;

figure;
imagesc(x,z,image_error);
colormap(gray);
set(gcf,'Position',[100 100 530 400]); 
set(gca,'Position',[.11 .03 .86 .85]);
set(gca,'xtick',0:0.2:2,'fontsize',13); 
set(gca,'ytick',0:0.2:2,'fontsize',13); 
set(gca,'linewidth',1.4);
xlabel('Distance (km)', 'fontsize',13);
ylabel('Depth (km)', 'fontsize',13);
set(gca,'xaxislocation','top');
caxis([-200 200]);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.3 4.0]);
print ./figures/figure11d.eps -depsc2 -r600;

%% Plot seismic traces located at X=0.3 km, X=0.5 km, X=0.7 km
X = 61;
trace_acoustic_x300 = image_ARTM_acoustic(:,X);
trace_viscoaco_x300 = image_ARTM_viscoac(:,X);
trace_q_comp_x300 = image_QRTM_viscoac(:,X);
X = 101;
trace_acoustic_x500 = image_ARTM_acoustic(:,X);
trace_viscoaco_x500 = image_ARTM_viscoac(:,X);
trace_q_comp_x500 = image_QRTM_viscoac(:,X);
X = 141;
trace_acoustic_x700 = image_ARTM_acoustic(:,X);
trace_viscoaco_x700 = image_ARTM_viscoac(:,X);
trace_q_comp_x700 = image_QRTM_viscoac(:,X);

figure;
subplot(3,1,1);
plot(z,trace_acoustic_x300,'r','LineWidth',1.5);
hold on;
plot(z,trace_viscoaco_x300,'b--','LineWidth',1.5);
hold on;
plot(z,trace_q_comp_x300,'g-.','LineWidth',1.5);
set(gcf,'Position',[100 100 600 500]); 
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',-4:4:4,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Depth (km)', 'fontsize',12);
xlim([0 1.1]);
ylim([-191 161]);
legend('Acoustic','Viscoacoustic','Compensated',...
      'Orientation','Horizontal','Location',[0.44 0.93 0.2 0.05]);
strings={'a)'};
annotation('textbox',[0.02,0.94,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);
subplot(3,1,2);
plot(z,trace_acoustic_x500,'r','LineWidth',1.5);
hold on;
plot(z,trace_viscoaco_x500,'b--','LineWidth',1.5);
hold on;
plot(z,trace_q_comp_x500,'g-.','LineWidth',1.5);
xlim([0 1.1]);
ylim([-201 301]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',-6:6:6,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Depth (km)', 'fontsize',12);
strings={'b)'};
annotation('textbox',[0.02,0.64,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);
subplot(3,1,3);
plot(z,trace_acoustic_x700,'r','LineWidth',1.5);
hold on;
plot(z,trace_viscoaco_x700,'b--','LineWidth',1.5);
hold on;
plot(z,trace_q_comp_x700,'g-.','LineWidth',1.5);
xlim([0 1.1]);
ylim([-181 181]);
set(gca,'xtick',0:0.2:2,'fontsize',12); 
set(gca,'ytick',-4:4:4,'fontsize',12); 
ylabel('Amplitude', 'fontsize',12);
xlabel('Depth (km)', 'fontsize',12);
strings={'c)'};
annotation('textbox',[0.01,0.34,0.2,0.05],'LineStyle','none',...
    'LineWidth',1,'String',strings,'fontsize',13);

% store the figure in the folder named 'figures'
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.0 5.0]);
print ./figures/figure12.eps -depsc2 -r600;
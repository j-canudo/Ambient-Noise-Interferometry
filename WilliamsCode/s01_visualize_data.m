clear all
close all

%datafile = './data/Fs1Hz_dx10m_Belgium_data.mat';
datafile = './data/Fs1Hz_dx10m_Tarifa_data.mat';
%datafile = './data/Fs1Hz_dx10m_A23_data.mat';

%Flag to save plots
save_plots = true;

%Parameters
nx = 1000; %number of channels
ns = 4200; %number of samples
fs = 1; %sampling rate (Hz)
dx = 10; %channel spacin (m)

%% Plot raw data
load(datafile);

%Create distance and time axis arrays
x = (0:nx-1)*dx;
t = (0:ns-1)/fs;

%Select a slice of the dataset
xmin = 0;
xmax = 500;
idx = (x>=xmin & x<=xmax);

tmin = 1000;
tmax = 1200;
idt = (t>=tmin & t<=tmax);

data_slice = all_data(idx,idt);
xx = x(idx);
tt = t(idt);

%Plot the data
load('RdBu.mat');
v=1e2;
fig1 = figure;
set(fig1, 'Position', [0 0 1200 800]); % set figure size
ax(1) = subplot(1,2,1);
pcolor(xx,tt,data_slice');
shading interp;
colormap(ax(1), RdBu);
clim([-v v]);
xlabel('Distance (m)');
ylabel('Time (s)');
title('Raw DAS data');

%% Plot frequency-wavenumber spectrum
%Take a larger slice of the dataset
xmin = 0;
xmax = 10000;
idx = (x>=xmin & x<=xmax);

tmin = 0;
tmax = 600;
idt = (t>=tmin & t<=tmax);

data_slice = all_data(idx,idt);
xx = x(idx);
tt = t(idt);

%Apply two FFTs to get the FK spectrum
fk = fft2(data_slice);
fk = fftshift(fk);
fk = 20*log10(abs(fk));

%Get FK axes
f = linspace(-fs/2,fs/2,length(tt));
k = linspace(-(1/dx)/2,(1/dx)/2,length(xx));

%Get theoretical dispersion relation
h = 20; %Water depth (m)
f_OSGW = sqrt(9.8*2*pi*k.*tanh(2*pi*k*h))*(0.5/pi);

%Plot FK spectrum
vmin = 60;
vmax = 120;
ax(2) = subplot(1,2,2);
pcolor(f,k,fk);
shading interp;
colormap(ax(2), 'Jet');
clim([vmin vmax]);
xlim([0 max(f)]);
hold on;
plot(f_OSGW,k,'k--');
xlabel('Distance (m)');
ylabel('Time (s)');
title('Raw DAS data');

% adjust spacing between subplots
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
pos = get(ax(1), 'Position');
pos(3) = pos(3) + 0.04;
set(ax(1), 'Position', pos);

%% Plot spectrum along cable
%For each channel, compute the spectrum
nns = floor(ns/2) + 1;
spec = zeros(nx,nns);
for ix=1:nx-1
    spec(ix,:) = 20*log10(abs(rfft(all_data(ix,:))));
end

%Get frequency axis
f = linspace(0,0.5,nns);

%Plot spectrum along cable
fig2 = figure;
set(fig2, 'Position', [0 0 1000 500]); % set figure size

ax = axes;
pcolor(ax, x, f, spec');
shading interp;
colormap(ax, 'jet');
clim(ax, [10 90]);
xlabel(ax, 'Distance (m)');
ylabel(ax, 'Frequency (Hz)');
title(ax, 'Spectrum vs. distance');

%% Save plots
if save_plots
    saveas(fig1,'./figs/01a_raw_data.png');
    saveas(fig2,'./figs/01b_spectrum.png');
end





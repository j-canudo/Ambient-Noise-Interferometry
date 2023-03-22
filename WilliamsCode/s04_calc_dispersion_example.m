clear all
close all

save_plots = true;

%Set path and source gather to plot
ncfdir = './ncfs/';
src = 100;

%Set parameters
nx = 200; %Number of channels in source gather
ns = 2^9; %Number of samples in source gather
fs = 1; %Sampling rate (Hz)
dx = 10; %Channel spacing (m)

x = (0:nx-1)*dx; %Offsets
t = linspace(-floor(ns/2),floor(ns/2)-1,ns)/fs; %Time lags

ncf_neg = zeros(nx,ns);
ncf_pos = zeros(nx,ns);

%Load and stack data
hrs = [0,1];
for hr_index=1:2
    hr = hrs(hr_index);
    fpath1 = strcat(ncfdir,'neg_ncf_src',num2str(src),'_hr',num2str(hr_index),'.mat');
    fpath2 = strcat(ncfdir,'pos_ncf_src',num2str(src),'_hr',num2str(hr_index),'.mat');
    load(fpath1);
    ncf_neg = ncf_neg + xc;
    load(fpath2);
    ncf_pos = ncf_pos + xc;
end

ncf_neg = ncf_neg / length(hrs);
ncf_pos = ncf_pos / length(hrs);

%Reorganize time lags
ncf_neg = cat(2,ncf_neg(:,floor(ns/2):end),ncf_neg(:,1:floor(ns/2)-1));
ncf_pos = cat(2,ncf_pos(:,floor(ns/2):end),ncf_pos(:,1:floor(ns/2)-1));
ncf_com = cat(2,ncf_neg(:,1:floor(ns/2)-1),ncf_pos(:,floor(ns/2):end));

%Parameters for dispersion image
vmin = 5; %Velocity range
vmax = 20;
vstep = 0.1;
fmin = 0.04; %Frequency range
fmax = 0.2;
xmax = 400; %Offset to calculate dispersion (m)
tmax = 128; %Time to calculate dispersion (s) (It also sets the frequency resolution)

%Cut positive lags
idx = (x<xmax);
idt = (t>=0 & t<=tmax);
trxc1 = ncf_com(idx,idt);

%Cut negative lags
idx = (x<xmax);
idt = (t<=0 & t>=-tmax);
trxc2 = ncf_com(idx,idt);
trxc2 = fliplr(trxc2);

%Compute dispersion images
[frq,vel,disp1] = calcDispersion2(trxc1,fs,dx,vmin,vmax,vstep,fmin,fmax);
[frq,vel,disp2] = calcDispersion2(trxc2,fs,dx,vmin,vmax,vstep,fmin,fmax);

%Pick dispersion curves
N = length(frq);
pick1 = zeros(1,N);
pick2 = zeros(1,N);
vals1 = zeros(1,N);
vals2 = zeros(1,N);

for i=1:N
    [vals(i),index1] = max(disp1(:,i));
    pick1(i) = vel(index1);
    [vals2(i), index2] = max(disp2(:,i));
    pick2(i) = vel(index2);
end

%Remove picks below a threshold
threshold = 2; %Manually set
bounds = (frq<fmax & frq>fmin); %Inside frequency range of interest
%Positive side
index1 = (vals1>threshold & bounds);
good_pick1 = pick1(index1);
good_freq1 = frq(index1);
bad_pick1 = pick1(~index1);
bad_freq1 = frq(~index1);
%Negative side
index2 = (vals2>threshold & bounds);
good_pick2 = pick2(index2);
good_freq2 = frq(index2);
bad_pick2 = pick2(~index2);
bad_freq2 = frq(~index2);


%Plot
load('RdBu.mat');
v = 1;
fig = figure();
set(fig, 'Position', [0, 0, 600, 1000]);
ax(1) = subplot(2,1,1);
pcolor(ax(1), -t(idt), x(idx), trxc1(:,end:-1:1));
shading(ax(1), 'flat');
colormap(ax(1), RdBu);
clim(ax(1), [-v, v]);
hold(ax(1), 'on');
pcolor(ax(1), t(idt), x(idx), trxc2(:,end:-1:1));
shading(ax(1), 'flat');
clim(ax(1), [-v, v]);
xlabel(ax(1), 'Time lag (s)');
ylabel(ax(1), 'Distance (m)');
xlim(ax(1), [-60, 60]);

ax(2) = subplot(2,1,2);
pcolor(ax(2), frq, vel, disp1);
shading(ax(2), 'flat');
colormap(ax(2), 'jet');
hold(ax(2), 'on');
pcolor(ax(2), -frq, vel, disp2);
shading(ax(2), 'flat');
xlabel(ax(2), 'Frequency (Hz)');
ylabel(ax(2), 'Phase speed (m/s)');
xlim(ax(2), [-0.3, 0.3]);
ylim(ax(2), [vmin, vmax]);

scatter(ax(2), good_freq1, good_pick1, 'k', 'filled', 'DisplayName', 'accepted');
scatter(ax(2), -good_freq2, good_pick2, 'k', 'filled');
scatter(ax(2), bad_freq1, bad_pick1, 'w', 'filled', 'DisplayName', 'rejected');
scatter(ax(2), -bad_freq2, bad_pick2, 'w', 'filled');
legend(ax(2));



%Save
if save_plots
    saveas(fig,'./figs/04_ex_disp.png');
end


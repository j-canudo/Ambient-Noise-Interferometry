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

%Plot data
v = 1;
load('RdBu.mat');
fig = figure;
set(fig, 'Position', [0 0 1500 600]); % set figure size

ax(1) = subplot(1,3,1);
pcolor(ax(1), x, t, ncf_neg');
shading interp;
colormap(ax(1), RdBu);
clim(ax(1), [-v v]);
xlabel(ax(1), 'Offset (m)');
ylabel(ax(1), 'Time lag (s)');
title(ax(1), 'Anti-causal (negative)');
xlim(ax(1), [0 1000]);
ylim(ax(1), [-100 100]);

ax(2) = subplot(1,3,2);
pcolor(ax(2), x, t, ncf_pos');
shading interp;
colormap(ax(2), RdBu);
clim(ax(2), [-v v]);
xlabel(ax(2), 'Offset (m)');
title(ax(2), 'Causal (positive)');
xlim(ax(2), [0 1000]);
ylim(ax(2), [-100 100]);

ax(3) = subplot(1,3,3);
pcolor(ax(3), x, t, ncf_com');
shading interp;
colormap(ax(3), RdBu);
clim(ax(3), [-v v]);
xlabel(ax(3), 'Offset (m)');
title(ax(3), 'Combined');
xlim(ax(3), [0 1000]);
ylim(ax(3), [-100 100]);

if save_plots
    saveas(fig,'./figs/03_ncf_stack.png');
end

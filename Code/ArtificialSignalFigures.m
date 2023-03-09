%% Initialization
clear all;
close all;
path = 'C:\Users\Usuario\Desktop\RawData\A23';
c = HDASdata(path);
c.getStrain_2D;
i=1;

%% Run after opening each file of same SNR
c.S_MS = S_MS_save;
c.c = c_save;
c.stacked_files = n_save;
c.calculatePWS;
SNR(i,:) = c.PWS(29,:);
i=i+1;


%% Plots for same SNR
for i=1:size(SNR,1)
    SNR_plot = SNR(i,:)/max(abs(SNR(i,:))) + 2*i;
    plot(linspace(-60,60,size(SNR,2)),SNR_plot);
    hold on;
end
y = [2 4 6 8 10 12 14];
yticks(y);
ylabels = {'1h','12h','24h','48h','7 days','15 days','31 days'};
yticklabels(ylabels);
ylabel('Stacked time');
xlabel('Offset time (s)');






%% ALL SNR REPRESENTATION
c.S_MS = S_MS_save;
c.c = c_save;
c.stacked_files = n_save;
c.calculatePWS;
SNR_general(i,:) = c.PWS(29,:);
i = i+1;

%% Plot for all SNR
for i=1:size(SNR_general,1)
    SNR_plot = SNR_general(i,:)/max(abs(SNR_general(i,:))) + 2*i;
    plot(linspace(-60,60,size(SNR_general,2)),SNR_plot);
    hold on;
end
y = [2 4 6 8 10 12 14];
yticks(y);
ylabels = {'-10','-5','-3','0','3','5','10'};
yticklabels(ylabels);
ylabel('SNR (dB)');
xlabel('Offset time (s)');



clear all;
close all;
tic;
path = 'C:\Users\Usuario\Desktop\RawData\A23';
c = HDASdata(path);
c.time_to_correlate = 1;
c.time_to_stack = 60;
total_files = 1440;
N_Channels = 2000;
PWS_matrix = zeros(N_Channels/2,c.time_to_correlate*250*2-1);

for i=100:10:N_Channels
    for j=1:total_files
        fprintf('Opening file %d of %d in iteration %d of %d\n',j,total_files/c.number_of_files,i/2,N_Channels/2);
        c.getStrain2D;
        c.choppData(i,i+60);
        c.highPassFilter(0.1,2000);
        c.temporalNormalization;
        c.smoothData(2);
        c.spectralWhitening;
        c.smoothData(2);
        c.keepStrainChannelsPair(size(c.Strain2D,1)/2,size(c.Strain2D,1)/2+6);
        c.correlateAndStack;
    end
    PWS_matrix(i/2,:) = c.PWS(2,:);
    save("A23PWS_Tomography_save.mat",'PWS_matrix');
    c.setFileConfiguration(1,1);
    c.resetPWS;
    c.resetStrain;
end
timer = toc;
fprintf('Total time: %f hours.\n',timer/3600);
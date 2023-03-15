tic
clear all;
close all;
path = 'D:\';
c = HDASdata(path);
c.setFileConfiguration(1,1);
c.time_to_correlate = 1;
c.time_to_stack = 60;
total_files = 16*60;

for i=1:total_files/c.number_of_files
    fprintf('Iteration %d of %d\n',i,total_files/c.number_of_files);
    c.getStrain2D;
    c.choppData(270,330);
    %c.downsampleData(5);
    c.highPassFilter(0.05,2000);
    %c.filterFK([],[],[2 50],'both');
    c.temporalNormalization;
    c.smoothData(2);
    c.spectralWhitening;
    c.smoothData(2);
    c.choppData(size(c.Strain2D,1)/2,size(c.Strain2D,1));
    c.correlateAndStack;
    c.findInfty;
    %c.savePWSdata(500);
end
c.cleanZeroRow;
c.plotPWS;
c.plotPWS_perChannel;

c.calculateDispersionVelocity([500 3000],[]);
c.plotDispersionVelocity;
c.displayComputationTime;
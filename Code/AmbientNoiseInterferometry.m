tic
clear all;
close all;
path = 'F:\Submarino\Tanda5';
c = HDASdata(path);
c.setFileConfiguration(1,45);
c.time_to_correlate = 60*45;
c.time_to_stack = c.time_to_correlate;
total_files = 45*2;

for i=1:total_files/c.number_of_files
    fprintf('Iteration %d of %d\n',i,total_files/c.number_of_files);
    c.getStrain2D;
    c.choppData(300,330);
    c.downsampleData(5);
    c.highPassFilter(0.1,[]);
    c.filterFK([],[],[2 50],'both');
%     c.temporalNormalization;
%     c.smoothData(2);
%     c.spectralWhitening;
%     c.smoothData(2);
%     c.choppData(size(c.Strain2D,1)/2,size(c.Strain2D,1));
    c.correlateAndStack;
    c.findInfty;
    c.savePWSdata(500);
end
%c.cleanZeroRow;
c.plotPWS;
c.plotPWS_perChannel;

c.calculateDispersionVelocity([5 20]);
c.plotDispersionVelocity;
c.displayComputationTime;
clear all;
%close all;
path = 'D:\';
save_path = 'C:\Users\Usuario\Desktop\JCanudo\Resultados\VariablesClasificar';
c = HDASdata(path);
c.setFileConfiguration(1,30);
c.time_to_correlate = 60*45;
c.time_to_stack = 60*45;
total_files = 9000;

for i=1:total_files/c.number_of_files
    fprintf('Iteration %d of %d\n',i,total_files/c.number_of_files);
    c.getStrain2D;
    c.choppData(300,500);
    c.downsampleData(1);
    c.detrendData;
    c.bandPassFilter([0.03,0.3],[]);
    %c.FKfilter(5,50,'neg',[]);
    c.whiten(0.03,0.3);
    %c.temporalNormalization;
%     c.smoothData(2);
%     c.spectralWhitening;
%     c.smoothData(2);
%     c.choppData(size(c.Strain2D,1)/2,size(c.Strain2D,1));
    c.correlateAndStackStandard;
    c.correlateAndStackPWS;
    %c.findInfty;
    %c.savePWSdata(1000,save_path);
end
%c.savePWSdata(1,save_path);
disp('Finished correlating.');
%c.cleanZeroRow;
c.plotNCF;
c.plotNCF_perChannel;
c.plotPWS;
c.plotPWS_perChannel;

c.calculateDispersionVelocity([5 20],60);
c.plotDispersionVelocity;
c.displayComputationTime;
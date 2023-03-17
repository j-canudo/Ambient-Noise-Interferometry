clear all;
%close all;
path = 'F:\Submarino\Tanda5';
save_path = 'C:\Users\Usuario\Desktop\JCanudo\Resultados\VariablesClasificar';
c = HDASdata(path);
c.setFileConfiguration(1,70);
c.time_to_correlate = 60*10;
c.time_to_stack = c.time_to_correlate;
total_files = 10;

for i=1:total_files/c.number_of_files
    fprintf('Iteration %d of %d\n',i,total_files/c.number_of_files);
    c.getStrain2D;
    c.choppData(300,330);
    c.downsampleData(1);
    %c.highPassFilter(0.05,[]);
    %c.lowPassFilter(30,[]);
    c.FKfilter(5,50,'neg',[]);
    %c.filterFK([],[],[5 20],'both','keep');
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

c.calculateDispersionVelocity([5 20],[]);
c.plotDispersionVelocity;
c.displayComputationTime;
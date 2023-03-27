clear all
close all
path = 'D:\';
c = HDASdata2(path);
c.setANIConfiguration(1,512,512*30,0.5,30,50,10000);
total_files = 60*10;

for i=1:total_files/c.number_of_files
    fprintf('Iteration %d/%d\n',i,total_files/c.number_of_files);
    c.getStrain2D;
    c.choppData(250,450);
    c.correlateAndStack;
end
disp('Finish correlating');
c.plotPWS;
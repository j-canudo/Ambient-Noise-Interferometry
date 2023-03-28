clear all
close all
path = 'D:\';
c = HDASdata2(path);
c.setANIConfiguration(80,512,512*30,10,20,100,5000);
total_files = 60*24*15;

for i=1:total_files/c.number_of_files
    fprintf('Iteration %d/%d\n',i,total_files/c.number_of_files);
    c.getStrain2D;
    c.choppData(250,450);
    c.correlateAndStack;
end
disp('Finish correlating');
c.plotPWS;
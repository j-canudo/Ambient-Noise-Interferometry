tic
clear all;
close all;
path = 'C:\Users\jcanu\Desktop\Universidad\TFM\RawData\Submarino';
c = HDASdata(path);
c.number_of_files = 10;
total_files =60;

for i=1:total_files/c.number_of_files
    fprintf('Iteration %d of %d\n',i,total_files/c.number_of_files);
    c.getStrain2D;
    c.choppData(300,330);
    c.highPassFilter(0.1,[]);
    c.calculateFK;
    c.stackFK;
    c.saveFK(100);
end
c.FK = c.stackedFK;
c.plotFK;
timer = toc;
fprintf('Total time: %f hours.\n',timer/3600);
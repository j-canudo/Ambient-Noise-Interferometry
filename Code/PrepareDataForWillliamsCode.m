tic
clear all
close all
path = 'D:\';
c = HDASdata(path);
c.setFileConfiguration(1,30);
total_time = 24; %In hours
all_data = zeros(2000,1);
for i=1:(total_time*60)/c.number_of_files
    c.getStrain2D;
    c.choppData(200,2200);
    c.downsampleData(5);
    all_data = [all_data c.Strain2D];
end
all_data(:,1) = [];
save('Fs5Hz_dx10m_A23_data.mat','all_data');
toc
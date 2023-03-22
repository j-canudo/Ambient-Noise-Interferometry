clear all
close all
tic
%datafile = './data/Fs1Hz_dx10m_Belgium_data.mat';
%datafile = './data/Fs1Hz_dx10m_Tarifa_data.mat';
datafile = './data/Fs1Hz_dx10m_A23_data.mat';
outdir = './ncfs/';

%Set parameters
fs = 1;
dx = 10;
fmin = 0.03;
fmax = 5;
nns = 2^9; %number of samples in a subwindow
nnw = nns/2+1; %number of samples in RFFT
nov = floor(nns/2);

%Set geometry
nx = 1000;
ddx = 10; %Spatial interval for sources
nr = 200;
src_list = 1:dx:nx-nr;

%Read data
load(datafile);
disp('Data loaded.');
toc

%Split data into two partially overlapping windows
ns = 2800; %Number of samples in each window
nwin = floor(size(all_data,2)/(ns/2))-1;
for i=1:nwin
    datasets(:,:,i) = all_data(:,(i-1)*(ns/2)+1:(i-1)*(ns/2)+ns);
end
nf = size(datasets,3);

%Set directions to filter
sgns = ["neg" "pos"];

%For each hour of data
for i=1:nf
    fprintf('Hour %d/%d.\n',i,nf);
    %Get data and dimensions
    data = datasets(:,:,i);
    nx = size(data,1);
    ns = size(data,2);
    nwn = floor(ns/nov)-1;

    %For each wavefield direction
    for sgn_index=1:length(sgns)
        sgn = sgns(sgn_index);
        fprintf('\t calculating NCFs for %s direction.\n',sgn);
        %Calculate cross-correlations for each source
        for src_index=1:length(src_list)
            src = src_list(src_index);
            rec_arr = src:src+nr-1;
            xc = get_src_gather_fk(data,src,rec_arr,nx,nns,nnw,nwn,nov,fmin,fmax,fs,dx,sgn);
            %Save
            filename = strcat(outdir,sgn,'_ncf_src',num2str(src-1),'_hr',num2str(i),'.mat');
            save(filename,"xc");
        end
    end
end

toc
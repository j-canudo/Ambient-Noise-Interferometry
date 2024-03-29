classdef HDASdata2 < handle
    properties
        data_path
        FileHeader
        Strain2D
        first_file = 1;
        number_of_files = 1;
        files_dir
        Spatial_Sampling_Meters = 0;
        N_Processed_Points = 0;
        Trigger_Frequency
        FiberRefStart = 0;
        FiberRefStop = 0;
        signal_type = 0; %1 if artificial signal
        is_chopped = 0; %1 if file is previously chopped
        StartFiberPoint = 0;
        EndFiberPoint = 0;
        N_Time_Samples = 0;
        chopped_in_program = 0;
        Data_Matrix
        Reference_Matrix
        FFT2D
        gamma = 1;
        PWS
        NCF
        NCF_phase
        stacked_files = 0;
        timeInStrain2D
        FK
        time_to_correlate = 60;
        time_to_stack = 60;
        processed_files = 0;
        stackedFK
        dispersionVelocity2D
        maxDispersionVelocity
        w
        c
        max_dispersion_velocities
        cmapHDASdata
        nSamplesCorrelate
        nSamplesStack
        fmin
        fmax
        vmin
        vmax
    end

    methods (Access = public)
        function obj = HDASdata2(path)
            tic;
            obj.data_path = path;
            obj.files_dir = dir(path);
            obj.files_dir = obj.files_dir(~ismember({obj.files_dir.name},{'.','..','$RECYCLE.BIN','System Volume Information'}));
            load('RdBu.mat','RdBu');
            obj.cmapHDASdata = RdBu;
        end

        function obj = getFileHeader(obj, fullPath)
            type = 'float64=>float64';
            fileID = fopen(fullPath);
            fHeaderSize = fread(fileID,1,type);
    
            header = zeros(fHeaderSize,1);
            header(1) = fHeaderSize;
            header(2:fHeaderSize) = fread(fileID, fHeaderSize-1,type);
            fclose(fileID);
            obj.FileHeader = header;
        end

        function obj =  getMeasurementSettings(obj,fullPath)
            obj.getFileHeader(fullPath);
            if obj.FileHeader(102) ~= 0
                print('ERROR: File is not "HDAS_2DRawData_Strain".');
                quit
            end

            if obj.FileHeader(175) == 1
                obj.is_chopped = 1;
            end

            if obj.FileHeader(176) == 1
                obj.signal_type = 1;
            end

            if obj.signal_type == 1
                obj.StartFiberPoint = 0;
                obj.EndFiberPoint = obj.FileHeader(2);
                obj.N_Processed_Points = obj.FileHeader(2);
                obj.FiberRefStart = 0;
                obj.FiberRefStop = 0;
                obj.Trigger_Frequency = obj.FileHeader(4);
                obj.Spatial_Sampling_Meters = obj.FileHeader(3);
            else
                obj.FiberRefStart = obj.FileHeader(18);
                obj.FiberRefStop = obj.FileHeader(20);
                obj.Spatial_Sampling_Meters = obj.FileHeader(2);
                obj.Trigger_Frequency = obj.FileHeader (7) / obj.FileHeader (16) / obj.FileHeader(99);
                if obj.is_chopped == 1
                    obj.StartFiberPoint = obj.FileHeader(13);
                    obj.EndFiberPoint = obj.FileHeader(15)-(obj.FiberRefStop-obj.FiberRefStart+1);
                else
                    if obj.EndFiberPoint == 0
                        obj.EndFiberPoint = obj.FileHeader(15);
                        obj.StartFiberPoint = obj.FileHeader(13);
                    end
                        obj.N_Processed_Points = obj.FileHeader(15) - obj.FileHeader(13);
                end
            end
        end

        function obj = choppData(obj, StartPoint, EndPoint)
            obj.StartFiberPoint = StartPoint;
            obj.EndFiberPoint = EndPoint-1;
            obj.chopped_in_program = 1;
            obj.Strain2D = obj.Strain2D(obj.StartFiberPoint:obj.EndFiberPoint,:);
            obj.N_Processed_Points = size(obj.Strain2D,1);
        end

        function obj = getStrainFromFile(obj, fullPath)
            fileID = fopen(fullPath);
            type = 'int16=>int16';
            raw_file_bad = fread(fileID, type);
            raw_file = raw_file_bad(801:end,1);
            fclose(fileID);
            
            mk_to_Strain = 1;
            if obj.signal_type
                obj.N_Time_Samples = obj.Trigger_Frequency*60;
            else
                obj.N_Time_Samples = length(raw_file)/(obj.FileHeader(15)-obj.FileHeader(13));
            end
            Traces_Matrix = zeros(obj.N_Processed_Points, obj.N_Time_Samples);
            for i=0:(obj.N_Time_Samples-1)
                trace_start = i*obj.N_Processed_Points;
                Traces_Matrix(:,i+1) = mk_to_Strain*raw_file(trace_start+1:trace_start+obj.N_Processed_Points);
            end
            
            if obj.signal_type
                obj.Reference_Matrix = zeros(obj.FiberRefStop-obj.FiberRefStart,obj.N_Time_Samples);
            else
                obj.Reference_Matrix = Traces_Matrix(obj.FiberRefStart+1:obj.FiberRefStop,:);
            end
            obj.Data_Matrix = Traces_Matrix;
        end

        function obj = getStrain2D(obj)
            fullPath = strcat(obj.data_path,'\',obj.files_dir(obj.first_file).name);
            obj.getMeasurementSettings(fullPath);
            
            Strain_2D_full = zeros(obj.EndFiberPoint-obj.StartFiberPoint,1);
            Strain_References = zeros(obj.FiberRefStop-obj.FiberRefStart+1,1);
            for i=1:obj.number_of_files
                disp(fullPath);
                obj.getStrainFromFile(fullPath);
                if i == 1
                    Strain_2D_full = zeros(size(obj.Data_Matrix,1),size(obj.Data_Matrix,2)*obj.number_of_files);
                    size_Strain_2D = size(obj.Data_Matrix,2);
                    Strain_References = zeros(size(obj.Reference_Matrix,1),size(obj.Reference_Matrix,2)*obj.number_of_files);
                    size_Strain_References = size(obj.Reference_Matrix,2);
                end
                Strain_2D_full(:,1+(i-1)*size_Strain_2D:(i-1)*size_Strain_2D+size_Strain_2D) = obj.Data_Matrix;
                Strain_References(:,1+(i-1)*size_Strain_References:(i-1)*size_Strain_References+size_Strain_References) = obj.Reference_Matrix;
                
                if i~=obj.number_of_files
                    fullPath = strcat(obj.data_path,'\',obj.files_dir(obj.first_file+i).name);
                end
            end
            obj.N_Processed_Points = size(Strain_2D_full,1);
            if obj.signal_type
                obj.Strain2D = Strain_2D_full/32767;
            else
                obj.Strain2D = obj.denoiseStrain_2D(Strain_2D_full,Strain_References);
            end
    
            obj.first_file = obj.first_file+obj.number_of_files;
            obj.timeInStrain2D = obj.N_Time_Samples/obj.Trigger_Frequency*obj.number_of_files;
        end

        function [Strain_2D_denoised] = denoiseStrain_2D(obj, Strain_Matrix, Strain_References)
            Sqrt_Ref_Fiber_Length = floor(sqrt(obj.FiberRefStop-obj.FiberRefStart+1));
            Segmented_Mean = zeros(1,Sqrt_Ref_Fiber_Length);

            Strain_Matrix_Variation=Strain_Matrix;
            Strain_Matrix_StrainReference=zeros(size(Strain_Matrix,1),1);
            Strain_Matrix_PreviousStrain=zeros(size(Strain_Matrix,1),1);
            Strain_Matrix_RefUpdates=zeros(size(Strain_Matrix,1),size(Strain_Matrix,2));

            Reference_Matrix_Variation=Strain_References;
            Reference_Matrix_StrainReference=zeros(size(Strain_References,1),1);
            Reference_Matrix_PreviousStrain=zeros(size(Strain_References,1),1);
            Reference_Matrix_RefUpdates=zeros(size(Strain_References,1),size(Strain_References,2));

            for i=1:size(Strain_Matrix,2) %For each temporal point
                %For Strain_2D
                for j=1:size(Strain_Matrix_Variation,1) % For each Spatial Point
                    if Strain_Matrix_Variation(j,i)>12000
                        Strain_Matrix_RefUpdates(j,i)=1; % Reference update flaged
                    end
                    Strain_Matrix_Variation(j,i)=(Strain_Matrix_Variation(j,i)-20000*Strain_Matrix_RefUpdates(j,i))/400; % Variation in digital samples with respect to the reference Fiber;
                    Strain_Matrix(j,i)=Strain_Matrix_StrainReference(j)+Strain_Matrix_Variation(j,i);
                    Strain_Matrix_StrainReference(j)=Strain_Matrix_StrainReference(j)+Strain_Matrix_Variation(j,i)*Strain_Matrix_RefUpdates(j,i);
                end
                %For Strain_References
                for j=1:size(Reference_Matrix_Variation,1) % For each Spatial Point
                    if Reference_Matrix_Variation(j,i)>12000
                        Reference_Matrix_RefUpdates(j,i)=1; % Reference update flaged
                    end
                    Reference_Matrix_Variation(j,i)=(Reference_Matrix_Variation(j,i)-20000*Reference_Matrix_RefUpdates(j,i))/400; % Variation in digital samples with respect to the reference fibre;
                    Strain_References(j,i)=Reference_Matrix_StrainReference(j)+Reference_Matrix_Variation(j,i);
                    Reference_Matrix_StrainReference(j)=Reference_Matrix_StrainReference(j)+Reference_Matrix_Variation(j,i)*Reference_Matrix_RefUpdates(j,i);
                end
                % Strain variation is defined as 0 for the first iteration; No laser phase noise or peak hopping is compensated in i=1;
                if i==1
                    Strain_Matrix_PreviousStrain=Strain_Matrix(:,1);
                    Reference_Matrix_PreviousStrain=Strain_References(:,1);
                end

                % Use Segmented Median of Mean
                for j=1:Sqrt_Ref_Fiber_Length
                    Segmented_Mean(j)= mean(Strain_References(1+(j-1)*Sqrt_Ref_Fiber_Length:(j)*Sqrt_Ref_Fiber_Length,i))-mean(Reference_Matrix_PreviousStrain(1+(j-1)*Sqrt_Ref_Fiber_Length:(j)*Sqrt_Ref_Fiber_Length));
                end

                Strain_Matrix(:,i)=Strain_Matrix(:,i)-median(Segmented_Mean);
                Strain_References(:,i)=Strain_References(:,i)-median(Segmented_Mean);

                Peak_Search_Window_Size=(obj.FileHeader(6)*8/(1+obj.FileHeader(10))); % Peak_Search_Window_Size = SampleRate*8/ChirpSlope;

                for j=1:size(Strain_Matrix_Variation,1)
                    if abs(Strain_Matrix(j,i)-Strain_Matrix_PreviousStrain(j))>Peak_Search_Window_Size/4
                        Strain_Matrix_RefUpdates(j,i)=1; % If discontinuity bigger than 1/4 Peak_Search_Window_Size after laser denoising; indicating peak hopping) => Reference update is flaged
                    end
                end

                Strain_Matrix_PreviousStrain=Strain_Matrix(:,i);
                Reference_Matrix_PreviousStrain=Strain_References(:,i);
            end
            Strain_2D_denoised = Strain_Matrix*10*2.88; %Final denoised strain matrix
        end

        function obj = plotStrain(obj, HDASpoint)
            if isempty(obj.Strain2D)
                obj.getStrain2D;
            end
            x_axis = linspace(0,size(obj.Strain2D,2)/obj.Trigger_Frequency,size(obj.Strain2D,2));
            figure(1);
            plot(x_axis,obj.Strain2D(HDASpoint,:));
            xlabel('Time (s)');
            ylabel('n\epsilon');
        end

        function obj = FFT(obj)
            obj.FFT2D = fftshift(fft(obj.Strain2D,[],2),2);
        end

        function obj = plotFFT(obj, HDASpoint)
            if isempty(obj.FFT2D)
                obj.FFT;
            end
            x_axis = linspace(-obj.Trigger_Frequency/2,obj.Trigger_Frequency/2,size(obj.Strain2D,2));
            FFT_abs = abs(obj.FFT2D(HDASpoint,:));
            figure(2);
            plot(x_axis,FFT_abs);
            xlabel('Frequency (Hz)');
            ylabel('Amplitude (n\epsilon/Hz)');
        end

        function obj = lowPassFilter(obj,fc,vc)
            if ~isempty(fc)
                [b,a] = butter(3,fc/(obj.Trigger_Frequency/2),'low');
                for i=1:size(obj.Strain2D,1)
                    obj.Strain2D(i,:) = filtfilt(b,a,obj.Strain2D(i,:));
                end
            end
            if ~isempty(vc)
                [b,a] = butter(3,(1/vc)/(1/obj.Spatial_Sampling_Meters/2),'low');
                for i=1:size(obj.Strain2D,2)
                    obj.Strain2D(:,i) = filtfilt(b,a,obj.Strain2D(:,i));
                end
            end
        end

        function obj = highPassFilter(obj,fc,vc)
            if ~isempty(fc)
                [b,a] = butter(3,fc/(obj.Trigger_Frequency/2),'high');
                for i=1:size(obj.Strain2D,1)
                    obj.Strain2D(i,:) = filtfilt(b,a,obj.Strain2D(i,:));
                end
            end
            if ~isempty(vc)
                [b,a] = butter(3,(1/vc)/(1/obj.Spatial_Sampling_Meters/2),'high');
                for i=1:size(obj.Strain2D,2)
                    obj.Strain2D(:,i) = filtfilt(b,a,obj.Strain2D(:,i));
                end
            end
        end

        function obj = bandPassFilter(obj,fc,vc)
            if ~isempty(fc)
                [b,a] = butter(3,[fc(1) fc(2)]/(obj.Trigger_Frequency/2),'bandpass');
                for i=1:size(obj.Strain2D,1)
                    obj.Strain2D(i,:) = filtfilt(b,a,obj.Strain2D(i,:));
                end
            end
            if ~isempty(vc)
                [b,a] = butter(3,([1/vc(1) 1/vc(2)])/(1/obj.Spatial_Sampling_Meters/2),'bandpass');
                for i=1:size(obj.Strain2D,2)
                    obj.Strain2D(:,i) = filtfilt(b,a,obj.Strain2D(:,i));
                end
            end
        end

        function obj = setFileConfiguration(obj, initial_file, file_number)
            obj.first_file = initial_file;
            obj.number_of_files = file_number;
        end

        function obj = plotStrainRows(obj,initial_HDAS_point, final_HDAS_point)
            x_axis = linspace(0,size(obj.Strain2D,2)/obj.Trigger_Frequency,size(obj.Strain2D,2));
            figure(3);
            for i=initial_HDAS_point:final_HDAS_point
                strain_data = obj.Strain2D(i,:)/max(abs(obj.Strain2D(i,:))) + i;
                plot(x_axis,strain_data);
                hold on;
            end
            hold off;
            ylim([0 final_HDAS_point+1]);
            xlabel('Time (s)');
            ylabel('HDAS point');
            title('Normalized n\epsilon');
        end

        function obj = temporalNormalization(obj)
            [env,~] = envelope(obj.Strain2D.');
            norm_env = env/size(env,1);
            obj.Strain2D = obj.Strain2D./(norm_env.');
            max_value = max(abs(obj.Strain2D),[],2);
            obj.Strain2D = bsxfun(@rdivide,obj.Strain2D,max_value);
            for i=1:obj.N_Processed_Points
                obj.Strain2D(i,:) = movmean(obj.Strain2D(i,:),5);
            end
            max_value = max(abs(obj.Strain2D),[],2);
            obj.Strain2D = bsxfun(@rdivide,obj.Strain2D,max_value);
%             mask = obj.Strain2D > 0;
%             obj.Strain2D(mask) = 1;
%             obj.Strain2D(~mask) = -1;
        end

        function obj = savePWSdata(obj,file_number,path)
            if mod(obj.processed_files,file_number) == 0
                PWS = obj.PWS;
                n = obj.processed_files;
                file_name = strcat(path,'\A23PWS_',num2str(n),'files.mat');
                save(file_name,'PWS','n');
            end
        end

        function obj = plotWaterfall(obj,waterfall_window_length,freq_min,freq_max)
            waterfall_matrix = zeros(obj.N_Processed_Points, obj.timeInStrain2D/waterfall_window_length);

            for j=1:obj.timeInStrain2D/waterfall_window_length
                FFT_matrix = abs(fft(obj.Strain2D(:,(j-1)*waterfall_window_length*obj.Trigger_Frequency+1:j*waterfall_window_length*obj.Trigger_Frequency),[],2));
                waterfall_matrix(:,j) = sum(FFT_matrix(:,freq_min*waterfall_window_length:freq_max*waterfall_window_length),2);
                clearvars FFT_matrix;
            end

            time_axis = linspace(0,obj.timeInStrain2D,obj.timeInStrain2D/waterfall_window_length);
            spatial_axis = linspace(obj.StartFiberPoint*obj.Spatial_Sampling_Meters,(obj.StartFiberPoint+obj.N_Processed_Points-1)*obj.Spatial_Sampling_Meters,obj.N_Processed_Points);

            figure(8);
            imagesc(time_axis,spatial_axis,waterfall_matrix);
            xlabel('Time (s)');
            ylabel('Fiber distance (m)');
            colorbar;
        end

        function obj = calculateFK(obj)
            %obj.FK = fftshift(fft2(obj.Strain2D,size(obj.Strain2D,1)*2,size(obj.Strain2D,2)*2));
            obj.FK = fftshift(fft2(obj.Strain2D));
            obj.FK = flip(obj.FK,1);
        end

        function obj = plotFK(obj)
            if isempty(obj.FK)
                obj.calculateFK;
            end
            k_axis = linspace(-1/obj.Spatial_Sampling_Meters/2,1/obj.Spatial_Sampling_Meters/2,obj.N_Processed_Points);
            f_axis = linspace(-obj.Trigger_Frequency/2,obj.Trigger_Frequency/2,obj.N_Time_Samples*obj.number_of_files);
            figure(9);
            imagesc(f_axis,k_axis,abs(obj.FK));
            xlabel('Frequency (Hz)');
            ylabel('Wave number (m^{-1})');
            title('FK');
            set(gca,'YDir','normal');
            clim([1e2 1e6]);
            colorbar;
        end

        function obj = findInfty(obj)
            x = sum(sum(isnan(obj.PWS)));
            if x ~= 0
                error('NaN found in PWS. Ending.');
            end
            x = sum(sum(isinf(obj.PWS)));
            if x ~= 0
                error('Inf found in PWS. Ending.');
            end
        end

        function obj = calculateDispersionVelocity(obj, vel, cut_time)
            if isempty(obj.PWS)
                error('No PWS calculated. Ending.');
            end
            if ~isempty(cut_time)
                NCF_cut = obj.PWS(:,(size(obj.PWS,2)+1)/2-cut_time*obj.Trigger_Frequency:(size(obj.PWS,2)+1)/2+cut_time*obj.Trigger_Frequency);
            else
                NCF_cut = obj.PWS;
            end

            obj.w = linspace(-obj.Trigger_Frequency/2,obj.Trigger_Frequency/2,size(NCF_cut,2));
            obj.c = linspace(vel(1),vel(2),1000);
            X_ij = fftshift(fft(NCF_cut,[],2),2);
            obj.dispersionVelocity2D = zeros(length(obj.w),length(obj.c));
            suma = 0;
            for i=1:length(obj.w)
                for j=1:length(obj.c)
                    for k=1:obj.N_Processed_Points
                        producto = (X_ij(k,i)/abs(X_ij(k,i)))*exp(-1i*2*pi*obj.w(i)*(k-1)*obj.Spatial_Sampling_Meters/obj.c(j));
                        suma = suma + producto;
                    end
                    obj.dispersionVelocity2D(i,j) = suma;
                    suma = 0;
                end
            end
            
            obj.max_dispersion_velocities = zeros(1,length(obj.w));
            for i=1:length(obj.w)
                [~,index] = max(abs(obj.dispersionVelocity2D(i,:)));
                obj.max_dispersion_velocities(1,i) = obj.c(index);
            end
        end

        function obj = plotDispersionVelocity(obj)
            figure(10);
            imagesc(obj.w,obj.c,abs(obj.dispersionVelocity2D.'));
            set(gca,'YDir','normal');
            xlim([obj.w(1) obj.w(end)]);
            ylim([obj.c(1) obj.c(end)]);
            xlabel('Frequency (Hz)');
            ylabel('Dispersion velocity (m/s)');
            title('Dispersion Velocity');
            hold on;
            plot(obj.w,obj.max_dispersion_velocities(1,:),'ro');
        end

        function obj = downsampleData(obj, new_frequency)
            if mod(obj.Trigger_Frequency,new_frequency) ~= 0
                x = round(obj.Trigger_Frequency/new_frequency);
                new_frequency = obj.Trigger_Frequency/x;
                fprintf('New Trigger Frequency set to %d in order to be a divisor of %d.\n',new_frequency, obj.Trigger_Frequency);
            end

            Strain2D_downsampled = zeros(size(obj.Strain2D,1),size(obj.Strain2D,2)/(obj.Trigger_Frequency/new_frequency));
            for i=1:obj.N_Processed_Points
                %Strain2D_downsampled(i,:) = downsample(obj.Strain2D(i,:),obj.Trigger_Frequency/new_frequency);
                Strain2D_downsampled(i,:) = obj.Strain2D(i,1:obj.Trigger_Frequency/new_frequency:end);
            end
            obj.Trigger_Frequency = new_frequency;
            obj.Strain2D = Strain2D_downsampled;
        end


        function obj = setANIConfiguration(obj, first_file, samplesToCorrelate, samplesToStack, fmin, fmax, vmin, vmax)
            obj.nSamplesCorrelate = samplesToCorrelate;
            obj.nSamplesStack = samplesToStack;
            obj.fmin = fmin;
            obj.fmax = fmax;
            obj.vmin = vmin;
            obj.vmax = vmax;
            obj.first_file = first_file;
            if isempty(obj.Trigger_Frequency)
                nfiles = ceil((samplesToStack/250)/60);
            else
                nfiles = ceil((samplesToStack/obj.Trigger_Frequency)/60);
            end
            obj.number_of_files = nfiles;
        end

        function obj = correlateAndStack(obj)
            obj.NCF = zeros(size(obj.Strain2D,1),obj.nSamplesCorrelate);
            obj.NCF_phase = zeros(size(obj.Strain2D,1),obj.nSamplesCorrelate);

            obj.stacked_files = obj.nSamplesStack/obj.nSamplesCorrelate;
            for i=1:obj.stacked_files
                xc = obj.getXCGather(obj.Strain2D(:,(i-1)*obj.nSamplesCorrelate+1:i*obj.nSamplesCorrelate),obj.fmin,obj.fmax,obj.vmin,obj.vmax,obj.Trigger_Frequency,obj.Spatial_Sampling_Meters ...
                    ,[]);
                obj.NCF = obj.NCF + xc;
                analytic = (hilbert(obj.NCF.')).';
                obj.NCF_phase = obj.NCF_phase + analytic./abs(analytic);
            end
            if isempty(obj.PWS)
                obj.PWS = zeros(size(obj.NCF));
            end
            stack = ((abs(obj.NCF_phase)/obj.stacked_files).^obj.gamma).*(obj.NCF/(obj.stacked_files));
            stack = cat(2,stack(:,floor(obj.nSamplesCorrelate/2):end),stack(:,1:floor(obj.nSamplesCorrelate/2)-1));
            obj.PWS = obj.PWS + stack;
            obj.first_file = obj.first_file-(obj.number_of_files-1);
        end

        function obj = plotPWS(obj)
            time_axis = linspace(-obj.nSamplesCorrelate/2/obj.Trigger_Frequency,obj.nSamplesCorrelate/2/obj.Trigger_Frequency,size(obj.PWS,2));
            offset_axis = linspace(0,obj.N_Processed_Points*obj.Spatial_Sampling_Meters,obj.N_Processed_Points);

            figure(5);
            imagesc(offset_axis,time_axis,obj.PWS.');
            colormap(obj.cmapHDASdata);
            colorbar;
            xlabel('Offset (m)');
            ylabel('Time (s)');
            title('PWS');
            set(gca,'YDir','normal');
            max_value = max(abs(obj.PWS(:)));
            clim([-max_value/5 max_value/5]);
            xlim([0 obj.N_Processed_Points*obj.Spatial_Sampling_Meters/2]);

            %Fk filter velocities display
            vmin = 1/obj.vmin*offset_axis;
            vmax = 1/obj.vmax*offset_axis;
            hold on;
            plot(offset_axis,vmin,'g');
            plot(offset_axis,vmax,'g');
        end

        function obj = plotSpectrum(obj)
            spec = zeros(size(obj.Strain2D)/2+1);
            for i=1:size(obj.Strain2D,1)
                spec(i,:) = 20*log10(abs(HDASdata2.rfft(obj.Strain2D(i,:))));
            end
            
            figure(4);
            set(gca, 'Position', [0 0 1000 500]);
            %Plot spectrum along cable
            ax = axes;
            pcolor(ax, linspace(0,size(obj.Strain2D,1)*obj.Spatial_Sampling_Meters,size(obj.Strain2D,1)), linspace(0,obj.Trigger_Frequency/2,size(obj.Strain2D,2)/2+1), spec');
            shading interp;
            colormap(ax, 'parula');
            ylim([0 20]);
            clim(ax, [10 90]);
            xlabel(ax, 'Distance (m)');
            ylabel(ax, 'Frequency (Hz)');
            title(ax, 'Spectrum vs. distance');
        end
    end

    methods(Static)
        function displayComputationTime
            timer = toc;
            hours = floor(timer/3600);
            timer = timer - hours*3600;
            minutes = floor(timer/60);
            seconds = round(timer - minutes*60);
            fprintf('Total computation time: %d hours %d minutes %d seconds.\n',hours,minutes,seconds);
        end

        function freqs = fftfreq(n,fs)
            freqs = linspace(-fs/2,fs/2,n);
        end

        function y = rfft(x)
            fft_x = fft(x);
            y = fft_x(1:(floor(length(fft_x)/2)+1));
        end

        function y = irfft(x)
            if mod(length(x),2) == 0
                even = true;
            else
                even = false;
            end
            if (even)
                n = 2 * (length(x) - 1 );
                s = length(x) - 1;
            else
                n = 2 * (length(x) - 1 );
                s = length(x);
            end
            xn = zeros(1,n);
            xn(1:length(x)) = x;
            xn(length(x):n) = conj(x(s:-1:2));
            y  = real(ifft(xn));
        end

        function out = detrendData(in)
            x = 1:size(in,2);
            out = zeros(size(in));
            for i=1:size(in,1)
                y = in(i,:);
                m = (sum(x.*y)-sum(x)*sum(y)/numel(x))/(sum(x.^2)-sum(x)^2/numel(x));
                b = mean(y) - m*mean(x);
                out(i,:) = y - (m*x + b);
            end
        end

        function xc = getXCGather(data,fmin,fmax,vmin,vmax,fs,dx,sgn)
            for i=1:size(data,1)
                data(i,:) = HDASdata2.detrendData(data(i,:));
                data(i,:) = HDASdata2.bp(data(i,:),fmin,fmax,fs);
            end
            data = HDASdata2.fk_filt(data,fs,dx,vmin,vmax,sgn);
            sp = zeros(size(data,1),size(data,2)/2+1);
            for i=1:size(data,1)
                sp(i,:) = HDASdata2.rfft(data(i,:));
                sp(i,:) = HDASdata2.whiten2(sp(i,:),fs,fmin,fmax);
            end
            spxc = zeros(size(data,1),size(data,2)/2+1);
            for i=1:size(data,1)
                spxc(i,:) = conj(sp(i,:)) .* sp(1,:);
            end
            xc = zeros(size(data));
            for i=1:size(data,1)
                xc(i,:) = HDASdata2.irfft(spxc(i,:));
            end
        end

        function y = bp(x,freqmin,freqmax,fs)
            [b,a] = butter(4,[freqmin freqmax]/(fs/2),'bandpass');
            y = filtfilt(b,a,x);
        end

        function tr_out = fk_filt(tr_in,fs,dx,vmin,vmax,sgn)
            if isempty(vmin)
                cmin = 50;
            else
                cmin = vmin;
            end
            if isempty(vmax)
                cmax = 5000;
            else
                cmax = vmax;
            end
            Nx = size(tr_in,1);
            Ns = size(tr_in,2);
            f0 = -fs/2:fs/Ns:fs/2-fs/Ns;
            k0 = -1/dx/2:1/dx/Nx:1/dx/2-1/dx/2/Nx;
            ft2 = fftshift(fft2(tr_in));
            [F,K] = meshgrid(f0,k0);
            C = F./K;
            filt = zeros(size(ft2));
            if isempty(sgn)
                filt(((C > cmin) & (C < cmax)) | ((C < -cmin) & (C > -cmax))) = 1;
            else
                if isequal(sgn,'pos')
                    filt((C > cmin) & (C < cmax)) = 1;
                else
                    filt((C < -cmin) & (C > -cmax)) = 1;
                end
            end
            filt = imgaussfilt(filt,3);
            ft2f = ft2.*filt;
            tr_out = real(ifft2(fftshift(ft2f)));
        end

        function spc = whiten2(spc,fs,low,high)
            N = (length(spc)-1)*2;
            i1 = ceil(low/(fs/N));
            i2 = ceil(high/(fs/N));
            spc(1:i1) = cos(linspace(pi/2,pi,i1)).^2 .* exp(1i*angle(spc(1:i1)));
            spc(i1:i2) = exp(1i*angle(spc(i1:i2)));
            spc(i2:end) = cos(linspace(pi,pi/2,length(spc)-i2+1)).^2 .* exp(1i*angle(spc(i2:end)));
        end

        function y = tempNormalization(x)
            y = zeros(size(x));
            mask = x > 0;
            y(mask) = 1;
            mask = x < 0;
            y(mask) = -1; 
        end
    end
end
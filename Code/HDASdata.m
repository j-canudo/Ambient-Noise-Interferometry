classdef HDASdata < handle
    properties
        data_path
        FileHeader
        Strain2D
        first_file = 1;
        number_of_files = 1;
        files_dir
        Spatial_Sampling_Meters = 0;
        N_Processed_Points = 0;
        Trigger_Frequency = 0;
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
        S_MS
        phase_MS
        gamma = 1;
        PWS
        NCF_total
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
        cmap_PWS = [0	0	1
0.00775193798449612	0.00775193798449612	1
0.0155038759689922	0.0155038759689922	1
0.0232558139534884	0.0232558139534884	1
0.0310077519379845	0.0310077519379845	1
0.0387596899224806	0.0387596899224806	1
0.0465116279069767	0.0465116279069767	1
0.0542635658914729	0.0542635658914729	1
0.0620155038759690	0.0620155038759690	1
0.0697674418604651	0.0697674418604651	1
0.0775193798449612	0.0775193798449612	1
0.0852713178294574	0.0852713178294574	1
0.0930232558139535	0.0930232558139535	1
0.100775193798450	0.100775193798450	1
0.108527131782946	0.108527131782946	1
0.116279069767442	0.116279069767442	1
0.124031007751938	0.124031007751938	1
0.131782945736434	0.131782945736434	1
0.139534883720930	0.139534883720930	1
0.147286821705426	0.147286821705426	1
0.155038759689922	0.155038759689922	1
0.162790697674419	0.162790697674419	1
0.170542635658915	0.170542635658915	1
0.178294573643411	0.178294573643411	1
0.186046511627907	0.186046511627907	1
0.193798449612403	0.193798449612403	1
0.201550387596899	0.201550387596899	1
0.209302325581395	0.209302325581395	1
0.217054263565891	0.217054263565891	1
0.224806201550388	0.224806201550388	1
0.232558139534884	0.232558139534884	1
0.240310077519380	0.240310077519380	1
0.248062015503876	0.248062015503876	1
0.255813953488372	0.255813953488372	1
0.263565891472868	0.263565891472868	1
0.271317829457364	0.271317829457364	1
0.279069767441860	0.279069767441860	1
0.286821705426357	0.286821705426357	1
0.294573643410853	0.294573643410853	1
0.302325581395349	0.302325581395349	1
0.310077519379845	0.310077519379845	1
0.317829457364341	0.317829457364341	1
0.325581395348837	0.325581395348837	1
0.333333333333333	0.333333333333333	1
0.341085271317829	0.341085271317829	1
0.348837209302326	0.348837209302326	1
0.356589147286822	0.356589147286822	1
0.364341085271318	0.364341085271318	1
0.372093023255814	0.372093023255814	1
0.379844961240310	0.379844961240310	1
0.387596899224806	0.387596899224806	1
0.395348837209302	0.395348837209302	1
0.403100775193798	0.403100775193798	1
0.410852713178295	0.410852713178295	1
0.418604651162791	0.418604651162791	1
0.426356589147287	0.426356589147287	1
0.434108527131783	0.434108527131783	1
0.441860465116279	0.441860465116279	1
0.449612403100775	0.449612403100775	1
0.457364341085271	0.457364341085271	1
0.465116279069767	0.465116279069767	1
0.472868217054264	0.472868217054264	1
0.480620155038760	0.480620155038760	1
0.488372093023256	0.488372093023256	1
0.496124031007752	0.496124031007752	1
0.503875968992248	0.503875968992248	1
0.511627906976744	0.511627906976744	1
0.519379844961240	0.519379844961240	1
0.527131782945737	0.527131782945737	1
0.534883720930233	0.534883720930233	1
0.542635658914729	0.542635658914729	1
0.550387596899225	0.550387596899225	1
0.558139534883721	0.558139534883721	1
0.565891472868217	0.565891472868217	1
0.573643410852713	0.573643410852713	1
0.581395348837209	0.581395348837209	1
0.589147286821706	0.589147286821706	1
0.596899224806202	0.596899224806202	1
0.604651162790698	0.604651162790698	1
0.612403100775194	0.612403100775194	1
0.620155038759690	0.620155038759690	1
0.627906976744186	0.627906976744186	1
0.635658914728682	0.635658914728682	1
0.643410852713178	0.643410852713178	1
0.651162790697675	0.651162790697675	1
0.658914728682171	0.658914728682171	1
0.666666666666667	0.666666666666667	1
0.674418604651163	0.674418604651163	1
0.682170542635659	0.682170542635659	1
0.689922480620155	0.689922480620155	1
0.697674418604651	0.697674418604651	1
0.705426356589147	0.705426356589147	1
0.713178294573644	0.713178294573644	1
0.720930232558140	0.720930232558140	1
0.728682170542636	0.728682170542636	1
0.736434108527132	0.736434108527132	1
0.744186046511628	0.744186046511628	1
0.751937984496124	0.751937984496124	1
0.759689922480620	0.759689922480620	1
0.767441860465116	0.767441860465116	1
0.775193798449612	0.775193798449612	1
0.782945736434109	0.782945736434109	1
0.790697674418605	0.790697674418605	1
0.798449612403101	0.798449612403101	1
0.806201550387597	0.806201550387597	1
0.813953488372093	0.813953488372093	1
0.821705426356589	0.821705426356589	1
0.829457364341085	0.829457364341085	1
0.837209302325581	0.837209302325581	1
0.844961240310078	0.844961240310078	1
0.852713178294574	0.852713178294574	1
0.860465116279070	0.860465116279070	1
0.868217054263566	0.868217054263566	1
0.875968992248062	0.875968992248062	1
0.883720930232558	0.883720930232558	1
0.891472868217054	0.891472868217054	1
0.899224806201550	0.899224806201550	1
0.906976744186047	0.906976744186047	1
0.914728682170543	0.914728682170543	1
0.922480620155039	0.922480620155039	1
0.930232558139535	0.930232558139535	1
0.937984496124031	0.937984496124031	1
0.945736434108527	0.945736434108527	1
0.953488372093023	0.953488372093023	1
0.961240310077519	0.961240310077519	1
0.968992248062015	0.968992248062015	1
0.976744186046512	0.976744186046512	1
0.984496124031008	0.984496124031008	1
0.992248062015504	0.992248062015504	1
1	1	1
1	0.992063492063492	0.992063492063492
1	0.984126984126984	0.984126984126984
1	0.976190476190476	0.976190476190476
1	0.968253968253968	0.968253968253968
1	0.960317460317460	0.960317460317460
1	0.952380952380952	0.952380952380952
1	0.944444444444444	0.944444444444444
1	0.936507936507937	0.936507936507937
1	0.928571428571429	0.928571428571429
1	0.920634920634921	0.920634920634921
1	0.912698412698413	0.912698412698413
1	0.904761904761905	0.904761904761905
1	0.896825396825397	0.896825396825397
1	0.888888888888889	0.888888888888889
1	0.880952380952381	0.880952380952381
1	0.873015873015873	0.873015873015873
1	0.865079365079365	0.865079365079365
1	0.857142857142857	0.857142857142857
1	0.849206349206349	0.849206349206349
1	0.841269841269841	0.841269841269841
1	0.833333333333333	0.833333333333333
1	0.825396825396825	0.825396825396825
1	0.817460317460317	0.817460317460317
1	0.809523809523810	0.809523809523810
1	0.801587301587302	0.801587301587302
1	0.793650793650794	0.793650793650794
1	0.785714285714286	0.785714285714286
1	0.777777777777778	0.777777777777778
1	0.769841269841270	0.769841269841270
1	0.761904761904762	0.761904761904762
1	0.753968253968254	0.753968253968254
1	0.746031746031746	0.746031746031746
1	0.738095238095238	0.738095238095238
1	0.730158730158730	0.730158730158730
1	0.722222222222222	0.722222222222222
1	0.714285714285714	0.714285714285714
1	0.706349206349206	0.706349206349206
1	0.698412698412698	0.698412698412698
1	0.690476190476191	0.690476190476191
1	0.682539682539683	0.682539682539683
1	0.674603174603175	0.674603174603175
1	0.666666666666667	0.666666666666667
1	0.658730158730159	0.658730158730159
1	0.650793650793651	0.650793650793651
1	0.642857142857143	0.642857142857143
1	0.634920634920635	0.634920634920635
1	0.626984126984127	0.626984126984127
1	0.619047619047619	0.619047619047619
1	0.611111111111111	0.611111111111111
1	0.603174603174603	0.603174603174603
1	0.595238095238095	0.595238095238095
1	0.587301587301587	0.587301587301587
1	0.579365079365079	0.579365079365079
1	0.571428571428571	0.571428571428571
1	0.563492063492064	0.563492063492064
1	0.555555555555556	0.555555555555556
1	0.547619047619048	0.547619047619048
1	0.539682539682540	0.539682539682540
1	0.531746031746032	0.531746031746032
1	0.523809523809524	0.523809523809524
1	0.515873015873016	0.515873015873016
1	0.507936507936508	0.507936507936508
1	0.500000000000000	0.500000000000000
1	0.492063492063492	0.492063492063492
1	0.484126984126984	0.484126984126984
1	0.476190476190476	0.476190476190476
1	0.468253968253968	0.468253968253968
1	0.460317460317460	0.460317460317460
1	0.452380952380952	0.452380952380952
1	0.444444444444444	0.444444444444444
1	0.436507936507937	0.436507936507937
1	0.428571428571429	0.428571428571429
1	0.420634920634921	0.420634920634921
1	0.412698412698413	0.412698412698413
1	0.404761904761905	0.404761904761905
1	0.396825396825397	0.396825396825397
1	0.388888888888889	0.388888888888889
1	0.380952380952381	0.380952380952381
1	0.373015873015873	0.373015873015873
1	0.365079365079365	0.365079365079365
1	0.357142857142857	0.357142857142857
1	0.349206349206349	0.349206349206349
1	0.341269841269841	0.341269841269841
1	0.333333333333333	0.333333333333333
1	0.325396825396825	0.325396825396825
1	0.317460317460317	0.317460317460317
1	0.309523809523810	0.309523809523810
1	0.301587301587302	0.301587301587302
1	0.293650793650794	0.293650793650794
1	0.285714285714286	0.285714285714286
1	0.277777777777778	0.277777777777778
1	0.269841269841270	0.269841269841270
1	0.261904761904762	0.261904761904762
1	0.253968253968254	0.253968253968254
1	0.246031746031746	0.246031746031746
1	0.238095238095238	0.238095238095238
1	0.230158730158730	0.230158730158730
1	0.222222222222222	0.222222222222222
1	0.214285714285714	0.214285714285714
1	0.206349206349206	0.206349206349206
1	0.198412698412698	0.198412698412698
1	0.190476190476190	0.190476190476190
1	0.182539682539683	0.182539682539683
1	0.174603174603175	0.174603174603175
1	0.166666666666667	0.166666666666667
1	0.158730158730159	0.158730158730159
1	0.150793650793651	0.150793650793651
1	0.142857142857143	0.142857142857143
1	0.134920634920635	0.134920634920635
1	0.126984126984127	0.126984126984127
1	0.119047619047619	0.119047619047619
1	0.111111111111111	0.111111111111111
1	0.103174603174603	0.103174603174603
1	0.0952380952380952	0.0952380952380952
1	0.0873015873015873	0.0873015873015873
1	0.0793650793650794	0.0793650793650794
1	0.0714285714285714	0.0714285714285714
1	0.0634920634920635	0.0634920634920635
1	0.0555555555555556	0.0555555555555556
1	0.0476190476190477	0.0476190476190477
1	0.0396825396825397	0.0396825396825397
1	0.0317460317460317	0.0317460317460317
1	0.0238095238095238	0.0238095238095238
1	0.0158730158730159	0.0158730158730159
1	0.00793650793650791	0.00793650793650791
1	0	0]
    end

    methods
        function obj = HDASdata(path)
            tic;
            obj.data_path = path;
            obj.files_dir = dir(path);
            obj.files_dir = obj.files_dir(~ismember({obj.files_dir.name},{'.','..','$RECYCLE.BIN','System Volume Information'}));
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

        function obj = removeMean(obj)
            strain_mean = mean(obj.Strain2D,2);
            obj.Strain2D = obj.Strain2D-strain_mean;
        end

        function obj = hannData(obj)
            H = hann(size(obj.Strain2D,2)).';
            obj.Strain2D = bsxfun(@times,obj.Strain2D,H);
        end

        function obj = smoothData(obj,percentage)
            %Horizontal smooth
            smooth_window = ones(1,size(obj.Strain2D,2));
            hann_window = hann(ceil(size(obj.Strain2D,2)*percentage/100*2));

            smooth_window(1,1:floor((percentage/100)*size(obj.Strain2D,2))) = hann_window(1:floor(end/2),1);
            smooth_window(1,floor(((100-percentage)/100)*size(obj.Strain2D,2))+1:end) = hann_window(floor(end/2)+1:end,1);

            obj.Strain2D = bsxfun(@times,obj.Strain2D,smooth_window);
            %Vertical smooth
            smooth_window = ones(size(obj.Strain2D,1),1);
            hann_window = hann(ceil(size(obj.Strain2D,1)*percentage/100*2));

            smooth_window(1:floor((percentage/100)*size(obj.Strain2D,1)),1) = hann_window(1:floor(end/2),1);
            smooth_window(floor(((100-percentage)/100)*size(obj.Strain2D,1))+1:end,1) = hann_window(floor(end/2)+1:end,1);

            obj.Strain2D = bsxfun(@times,obj.Strain2D,smooth_window);
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

        function obj = spectralWhitening(obj)
            obj.FFT;
            smoothed_amplitudes = smoothdata(abs(obj.FFT2D),1,'movmean',20);
            obj.FFT2D = obj.FFT2D./smoothed_amplitudes;
            obj.FFT2D(isnan(obj.FFT2D)) = 0;
            obj.Strain2D = real(ifft(ifftshift(obj.FFT2D,2),[],2));
            max_value = max(abs(obj.Strain2D),[],2);
            obj.Strain2D = bsxfun(@rdivide,obj.Strain2D,max_value);
%             obj.FFT;
%             obj.FFT2D = obj.FFT2D./abs(obj.FFT2D);
%             obj.FFT2D(isnan(obj.FFT2D)) = 0;
%             obj.Strain2D = real(ifft(ifftshift(obj.FFT2D,2),[],2));
%             max_value = max(abs(obj.Strain2D),[],2);
%             obj.Strain2D = bsxfun(@rdivide,obj.Strain2D,max_value);
        end

        function obj = whiten(obj,fmin,fmax)
            obj.FFT;
            %sp = obj.rfft(obj.Strain2D);
            %N = (size(sp,2)-1)*2;
            N = size(obj.FFT2D,2);
            i1 = N/2+ceil(fmin/(obj.Trigger_Frequency/N));
            i2 = N/2+ceil(fmax/(obj.Trigger_Frequency/N));

            obj.FFT2D(:,N/2+1:i1) = cos(linspace(pi/2,pi,i1-N/2)).^2 .* exp(1i*angle(obj.FFT2D(:,N/2+1:i1)));
            obj.FFT2D(:,N-i1+1:N/2) = cos(linspace(pi,pi/2,i1-N/2)).^2 .* exp(1i*angle(obj.FFT2D(:,N-i1+1:N/2)));
            obj.FFT2D(:,i1:i2) = exp(1i*angle(obj.FFT2D(:,i1:i2)));
            obj.FFT2D(:,N-i2:N-i1) = exp(1i*angle(obj.FFT2D(:,N-i2:N-i1)));
            obj.FFT2D(:,i2:end) = cos(linspace(pi,pi/2,size(obj.FFT2D,2)-i2+1)).^2 .* exp(1i*angle(obj.FFT2D(:,i2:end)));
            obj.FFT2D(:,1:N-i2+1) = cos(linspace(pi/2,pi,size(obj.FFT2D,2)-i2+1)).^2 .* exp(1i*angle(obj.FFT2D(:,1:N-i2+1)));
            obj.Strain2D = real(ifft(ifftshift(obj.FFT2D,2),[],2));
        end

        function NCF = crossCorrelationTimeDomain(obj,strain)
            NCF = zeros(obj.N_Processed_Points,size(strain,2)*2-1);
            for i=1:obj.N_Processed_Points
                NCF(i,:) = xcorr(strain(1,:),strain(i,:));
            end
        end

        function NCF = crossCorrelationFrequencyDomain(obj,strain)
            NCF = zeros(size(strain,1),size(strain,2)*2-1);
            a = [zeros(1,length(strain(1,:))-1) strain(1,:)];
            FFT_CH1 = fft(a);
            for i=1:obj.N_Processed_Points
                b = [strain(i,:) zeros(1,length(strain(i,:))-1)];
                FFT_CH2 = fft(b);
                NCF(i,:) = real(ifft(FFT_CH1.*conj(FFT_CH2)));
            end
        end

        function obj = stackPWS(obj,NCF)
            if isempty(obj.S_MS)
                obj.S_MS = zeros(size(NCF));
            end
            if isempty(obj.phase_MS)
                obj.phase_MS = zeros(size(NCF));
            end
            analytic = (hilbert(NCF.')).';
            phase = analytic./abs(analytic);
            %phase = angle(analytic);
            %phase = unwrap(phase);

            obj.S_MS = obj.S_MS + NCF;
            %obj.phase_MS = obj.phase_MS + exp(1i*phase);
            obj.phase_MS = obj.phase_MS + phase;
        end

        function obj = correlateAndStackStandard(obj)
            NCF = zeros(obj.N_Processed_Points,size(obj.Strain2D,2)*2-1);
            if isempty(obj.PWS)
                obj.NCF_total = zeros(size(NCF));
            end
            for i=1:obj.N_Processed_Points
                NCF(i,:) = xcorr(obj.Strain2D(1,:),obj.Strain2D(i,:));
            end
            obj.NCF_total = obj.NCF_total + flip(NCF,2);
        end

        function obj = correlateAndStackPWS(obj)
            for j=1:obj.number_of_files/(obj.time_to_stack/60)
                for i=1:obj.time_to_stack/obj.time_to_correlate
                    strain = obj.Strain2D(:,(j-1)*obj.time_to_stack*obj.Trigger_Frequency+(i-1)*obj.Trigger_Frequency*obj.time_to_correlate+1:i*obj.Trigger_Frequency*obj.time_to_correlate+(j-1)*obj.time_to_stack*obj.Trigger_Frequency);
                    NCF = obj.crossCorrelationTimeDomain(strain);
                    %NCF = obj.crossCorrelationFrequencyDomain(strain);
                    obj.stackPWS(NCF);
                end
                obj.stacked_files = size(obj.Strain2D,2)/obj.N_Time_Samples;
                obj.processed_files = obj.processed_files + obj.stacked_files;
                obj.calculatePWS;
            end
        end

        function obj = calculatePWS(obj)
            if isempty(obj.PWS)
                obj.PWS = zeros(size(obj.S_MS));
            end
            x = ((abs(obj.phase_MS)/obj.stacked_files).^obj.gamma).*(obj.S_MS/(obj.stacked_files));
            obj.PWS = obj.PWS + flip(x,2);
            obj.phase_MS = [];
            obj.S_MS = [];
            obj.stacked_files = 0;
        end

        function obj = resetPWS(obj)
            obj.PWS = [];
            obj.S_MS = [];
            obj.phase_MS = [];
            obj.stacked_files = 0;
            obj.processed_files = 0;
        end

        function obj = savePWSdata_artificial(obj)
            if ismember(obj.processed_files,[60 720 1440 2880 10080 21600 44640])
                S_MS_save = obj.S_MS;
                c_save = obj.phase_MS;
                n_save = obj.processed_files;
                file_name = strcat('PWSdata_',num2str(n_save),'.mat');
                save(file_name,'S_MS_save','c_save','n_save');
            end
        end

        function obj = savePWSdata(obj,file_number,path)
            if mod(obj.processed_files,file_number) == 0
                PWS = obj.PWS;
                n = obj.processed_files;
                file_name = strcat(path,'\A23PWS_',num2str(n),'files.mat');
                save(file_name,'PWS','n');
            end
        end

        function obj = plotNCF(obj)
            time_axis = linspace(-obj.time_to_correlate,obj.time_to_correlate,size(obj.NCF_total,2));
            offset_axis = linspace(0,obj.N_Processed_Points*obj.Spatial_Sampling_Meters,obj.N_Processed_Points);

            figure(4);
            imagesc(offset_axis,time_axis,obj.NCF_total.');
            colormap(obj.cmap_PWS);
            colorbar;
            xlabel('Offset (m)');
            ylabel('Time (s)');
            title('NCF Standard');
            set(gca,'YDir','normal');
            clim([-1e3 1e3]);
        end

        function obj = plotPWS(obj)
            time_axis = linspace(-obj.time_to_correlate,obj.time_to_correlate,size(obj.PWS,2));
            offset_axis = linspace(0,obj.N_Processed_Points*obj.Spatial_Sampling_Meters,obj.N_Processed_Points);

            figure(5);
            imagesc(offset_axis,time_axis,obj.PWS.');
            colormap(obj.cmap_PWS);
            colorbar;
            xlabel('Offset (m)');
            ylabel('Time (s)');
            title('PWS');
            set(gca,'YDir','normal');
            clim([-1e3 1e3]);
        end

        function obj = plotNCF_perChannel(obj)
            time_axis = linspace(-obj.time_to_correlate,obj.time_to_correlate,size(obj.NCF_total,2));
            k=0;
            NCF_plot = zeros(size(obj.NCF_total));
            for i=1:size(obj.NCF_total,1)
                NCF_plot(i,:) = obj.NCF_total(i,:)/max(abs(obj.NCF_total(i,:)))*5+k*obj.Spatial_Sampling_Meters;
                k=k+1;
            end
            figure(6);
            for i=1:size(obj.NCF_total,1)
                plot(time_axis,NCF_plot(i,:));
                hold on;
            end
            hold off;
            xlabel('Time lag(s)');
            ylabel('Offset (m)');
            title('NCF Standard');
            set(gca,'YDir','normal');
            xlim([-obj.time_to_correlate obj.time_to_correlate]);
            ylim([-5 obj.N_Processed_Points*obj.Spatial_Sampling_Meters+5]);
        end

        function obj = plotPWS_perChannel(obj)
            time_axis = linspace(-obj.time_to_correlate,obj.time_to_correlate,size(obj.PWS,2));
            k=0;
            PWS_plot = zeros(size(obj.PWS));
            for i=1:size(obj.PWS,1)
                PWS_plot(i,:) = obj.PWS(i,:)/max(abs(obj.PWS(i,:)))*5+k*obj.Spatial_Sampling_Meters;
                k=k+1;
            end
            figure(7);
            for i=1:size(obj.PWS,1)
                plot(time_axis,PWS_plot(i,:));
                hold on;
            end
            hold off;
            xlabel('Time lag(s)');
            ylabel('Offset (m)');
            title('PWS');
            set(gca,'YDir','normal');
            xlim([-obj.time_to_correlate obj.time_to_correlate]);
            ylim([-5 obj.N_Processed_Points*obj.Spatial_Sampling_Meters+5]);
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

        function obj = FKfilter(obj,cmin,cmax,direction,keep)
            Nx = size(obj.Strain2D,1);
            Ns = size(obj.Strain2D,2);
            f0 = obj.fftfreq(Ns,obj.Trigger_Frequency);
            k0 = obj.fftfreq(Nx,1/obj.Spatial_Sampling_Meters);

            obj.calculateFK;
            [F,K] = meshgrid(f0,k0);
            C = F./K;
            filt = zeros(size(obj.FK));
            
            if isequal(direction,'pos')
                filt((C > cmin) & (C < cmax)) = 1;
            end
            if isequal(direction,'neg')
                filt((C < -cmin) & (C > -cmax)) = 1;
            end
            if isempty(direction)
                filt(((C > cmin) & (C < cmax)) || ((C < -cmin) & (C > -cmax))) = 1;
            end

            if isequal(keep,'discard')
                filt = double(~logical(filt));
            end

            filt = imgaussfilt(filt,3);
            obj.FK = obj.FK.*filt;
            obj.Strain2D = real(ifft2(ifftshift(obj.FK)));
        end

        function obj = filterFK(obj,freq,k,v,direction,keep)
            obj.calculateFK;
            mask = logical(ones(size(obj.FK)));
            k_axis = linspace(-1/obj.Spatial_Sampling_Meters/2,1/obj.Spatial_Sampling_Meters/2,size(obj.FK,1));
            f_axis = linspace(-obj.Trigger_Frequency/2,obj.Trigger_Frequency/2,size(obj.FK,2));
            v_matrix = k_axis.'*(1./f_axis);
            
            if ~isempty(v)
                v_min = 1/v(1);
                v_max = 1/v(2);

                mask1 = v_matrix<v_min;
                mask1 = mask1 & v_matrix>-v_min;
                mask2 = v_matrix>v_max;
                mask2 = mask2 | v_matrix<-v_max;
                mask = mask1 & mask2;
            end

            if ~isempty(freq)
                freq_min = freq(1);
                freq_max = freq(2);
                delta_f = obj.Trigger_Frequency/size(obj.FK,2);
                default_frequency_highpass = 0.04;

                mask(:,length(f_axis)/2-round(default_frequency_highpass/delta_f):length(f_axis)/2+round(default_frequency_highpass/delta_f)) = false; %By default, highpass filter >0.04Hz
                mask(:,length(f_axis)/2+freq_max/delta_f:end) = false; %Eliminate f>freq_max
                mask(:,1:length(f_axis)/2-freq_max/delta_f) = false;
                mask(:,length(f_axis)-freq_min/delta_f:length(f_axis)/2+freq_min/delta_f) = false; %Eliminate f<f_min
            end

            if ~isempty(k)
                k_min = k(1);
                k_max = k(2);
                delta_k = (1/obj.Spatial_Sampling_Meters)/size(obj.FK,1);

                mask(1:length(k_axis)/2-k_max/delta_k,:) = false; %Eliminate k>k_max
                mask(length(k_axis)/2+k_max/delta_k:end,:) = false;
                mask(length(k_axis)/2-k_min/delta_k:length(k_axis)/2+k_min/delta_k,:) = false; %Eliminate k<k_min
            end

            if isequal(direction,'normal')
                a = ones(size(obj.FK)/2);
                b = zeros(size(obj.FK)/2);
                mask = [b a;a b];
            end
            if isequal(direction,'reverse')
                a = ones(size(obj.FK)/2);
                b = zeros(size(obj.FK)/2);
                mask = [a b; b a];
            end
            mask = double(mask);

            filt = imgaussfilt(mask,3);
            mask = mask.*filt;
            if isequal(keep,'keep')
                obj.FK(~mask) = 0;
            end
            if isequal(keep,'discard')
                obj.FK(mask) = 0;
            end
            if ~isequal(keep,'keep') && ~isequal(keep,'discard')
                error('Last argument should be keep or discard.');
            end
            obj.FK = flip(obj.FK,1);
            recovered_strain = real(ifft2(ifftshift(obj.FK),size(obj.Strain2D,1),size(obj.Strain2D,2)));
            obj.Strain2D = recovered_strain(1:size(obj.Strain2D,1),1:size(obj.Strain2D,2));
        end

        function obj = normalizeFK(obj,freq,k)
            if isempty(obj.FK)
                obj.calculateFK;
            end
            if freq
                max_values = max(abs(obj.FK),[],1);
                max_values(1:length(max_values)/2) = 1;
                obj.FK = bsxfun(@rdivide,obj.FK,max_values);
            end
            if k
                max_values = max(abs(obj.FK),[],2);
                obj.FK = bsxfun(@rdivide,obj.FK,max_values);
            end
        end

        function obj = stackFK(obj)
            if isempty(obj.stackedFK)
                obj.stackedFK = zeros(size(obj.FK));
            end
            obj.stackedFK = obj.stackedFK + obj.FK;
            obj.stacked_files = obj.stacked_files + obj.number_of_files;
        end

        function obj = saveFK(obj,n_files)
            if mod(obj.stacked_files,n_files) == 0
                stackFK = obj.stackedFK;
                n = obj.stacked_files;
                save('stackFK.mat','stackFK','n');
            end
        end

        function obj = keepStrainChannelsPair(obj, first_channel, second_channel)
            x1 = obj.Strain2D(first_channel,:);
            x2 = obj.Strain2D(second_channel,:);
            obj.Strain2D = zeros(2,length(x1));
            obj.Strain2D(1,:) = x1;
            obj.Strain2D(2,:) = x2;
            obj.N_Processed_Points = 2;
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

        function obj = cleanZeroRow(obj)
            obj.PWS(:,round(size(obj.PWS,2)/2)) = zeros(size(obj.PWS,1),1);
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

        function obj = calculateDispersionVelocityFK(obj,w,vel)
            if isempty(obj.FK)
                error('No FK calculated. Ending.');
            end

            k_axis = linspace(-1/obj.Spatial_Sampling_Meters/2,1/obj.Spatial_Sampling_Meters/2,size(obj.FK,1));
            f_axis = linspace(-obj.Trigger_Frequency/2,obj.Trigger_Frequency/2,size(obj.FK,2));
            v_matrix = k_axis.'*(1./f_axis);
            delta_f = obj.Trigger_Frequency/size(obj.FK,2);
            delta_k = (1/obj.Spatial_Sampling_Meters)/size(obj.FK,1);
            delta_v = delta_f/delta_k;

            obj.c = [vel(1):delta_v:vel(2)];
            obj.w = [-w:delta_f:w];
            obj.dispersionVelocity2D = zeros(length(obj.c),length(obj.w));

            for i=1:length(obj.c)
                v_min = 1/obj.c(i);
                v_max = v_min + delta_v;

                mask1 = v_matrix<v_min;
                mask1 = mask1 & v_matrix>-v_min;
                mask2 = v_matrix>v_max;
                mask2 = mask2 | v_matrix<-v_max;
                mask = mask1 & mask2;

                x = obj.FK;
                x(~mask) = 0;

                for j=1:length(obj.w)
                    obj.dispersionVelocity2D(i,j) = sum(abs(x(:,size(x,2)/2+(obj.w(j)/delta_f))));
                end
            end

            obj.dispersionVelocity2D = obj.dispersionVelocity2D.';

            obj.max_dispersion_velocities = zeros(1,length(obj.w));
            for i=1:length(obj.w)
                [~,index] = max(abs(obj.dispersionVelocity2D(i,:)));
                obj.max_dispersion_velocities(1,i) = obj.c(index);
            end
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
            N = size(x,2);
            y = fft(x,[],2);
            if mod(N,2) == 0
                y = y(:,1:N/2+1);
            else
                y = y(:,1:(N+1)/2);
            end
        end

        function y = irfft(x)
            x_pad = zeros(size(x,1),size(x,2)-2);
            full_x = [x x_pad];
            y = ifft(full_x);
            y = real(y);
        end
    end
end
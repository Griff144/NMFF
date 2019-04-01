tic
clear all
% close all

PerformAutocorrelation = 0;
Reducing = 0;
datatype = 1; %0 for .bin and 1 for xls
plot_irf = 0;
plot_events = 0;

threshold = 0.45;%4.8;%0.8;%0.0055;
step = 0.01;
for i = 1:1000
    thresholds(i,1) = step;
    step = step + 0.01;
end
%% Load data
if datatype == 0
    Filename = '\\10.164.16.153\Data\Detrick\08-10-2018\10-Aug-2018_11_51_14Channel_B_Raw_Data.bin';%'\\10.164.16.153\Data\Scott  Griffin\deconvolution\DoxycyclineHyclate\08-Aug-2017_16_12_28Channel_B_Raw_Data.bin';
    fid = fopen(Filename);
    filter = 50;
    Range = 400; %input range for alazar card
    Buffers = 36; %amount of concatenated hits
    temp2 = zeros(2050000,Buffers);
    reduction = 1;
    Data = fread(fid,'uint16');%dlmread(Filename,'\t');
    %% Convert into voltage and make positive
    Rescaled_Data = ((Data / 2^15) * -Range) + Range;%converts to voltage, flips to positive, and shifts baseline towards 0
    split_data = reshape(Rescaled_Data,[],Buffers);
    for i = 1:Buffers;
     baseline2 = mean(split_data(1:40000,i)); %finds the mean of the first millisecond of data in each trace
     temp2(:,i) = (split_data(:,i)-baseline2); %subtracts the baseline mean from the corresponding trace and saves them
    end
%     baseline = mean(Rescaled_Data(1:40000/reduction));
%     Rescaled_Data = Rescaled_Data-baseline;
end

if datatype == 1
    Filename = '\\10.164.16.153\Data\Scott  Griffin\deconvolution\ROC plot\TimeTraces&GroundTruth\TimeTrace9.csv';
    Rescaled_Data = csvread(Filename);
    Filename2 = '\\10.164.16.153\Data\Scott  Griffin\deconvolution\ROC plot\true_irf.csv';
    GroundTruthIRF = csvread(Filename2,0,0);
    %Filename3 = 'C:\Users\scott\OneDrive\Documents\triboluminescence\deconvolution\PeakShape.csv';
    %PeakShape = csvread(Filename3,0,0);
    filter = 50;%100;
    reduction = 1;
    Buffers = 1;
    temp2 = Rescaled_Data;
end
%irf_struct = load('\\10.164.16.153\Data\CaseySmith\Deconvolution\IRFs.mat');
IRFs = xlsread('\\10.164.16.153\Data\Scott  Griffin\deconvolution\Gaus_IRF.xlsx');%csvread('\\10.164.16.153\Data\Scott  Griffin\deconvolution\IRF_result_matrix_random.csv')';%irf_struct([1]).outputIrf;
% IRFs = csvread('\\10.164.16.153\Data\Scott  Griffin\deconvolution\test_irf.csv');
IRF_test_results = csvread('\\10.164.16.153\Data\Scott  Griffin\deconvolution\IRF_result_matrix4.csv');

recovered_events = cell(1,Buffers);
recovered_irf = cell(1,Buffers);
for n = 1:Buffers
%     if value1(n) == 0
%         n = n + 1;
%     end
    Reduced_Data = temp2;%temp2(value1(n)-1000:value2(n)+1000,n);
%% Create baseline correcting high-pass filter and apply it
    x_hp = linspace(-3,3);
    pdf = normpdf(x_hp,0,1);
    pdf_normalized = -pdf/sum(pdf);
    impulse = zeros(1,100);
    for i = 1:length(x_hp)
        if i == 50
            impulse(1,i) = abs(sum(pdf_normalized));
        else
            impulse(1,i) = 0;
        end
    end
    final_filter = pdf_normalized + impulse;
    Reduced_Data = conv(final_filter,Reduced_Data);
    Reduced_Data = Reduced_Data(50:length(Reduced_Data)-50);
    saved_variable = Reduced_Data;
    %% Fit 
    filtersize = filter/reduction; %size of data transient 
    irf_results = zeros(filtersize,1);
    delta = 0;
    how_many = 10;
    x = linspace(.1,5,filtersize);
    %data_results = zeros(length(Reduced_Data),how_many);
    MuGuess = 1; %where the guess starts

    for k = 1:how_many
        MuGuess = MuGuess + delta; %original guess value for the eponential pdf
        comb = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]';%1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]';
        decay_guess = comb .* exppdf(x,MuGuess)'; %original guess pdf 
    %     decay_guess = exppdf(x,MuGuess)';
        keepgoing = 1; 
        condition1 = 1; %condition for X2(kinda) to be getting smaller
        condition2 = 1; %condition for the amount of recovered photon events to be getting smaller
        counter = 0;%for how many times it needs to iterate between recovering the irf and amplitudes
        counter2 = 0;
        tstart = tic;
        while keepgoing == 1 
        %% making the orignial E matrix
            data = Reduced_Data;
            L = filtersize;
            shift = 0;
            E_init = zeros(filtersize,filtersize); %change filtersize to L to make it like original
            for c = 1:L
                for r = 1:length(decay_guess)
                    if r+shift<L+1
                        E_init(r+shift,c) = decay_guess(r);
                    else
                        E_init(r+shift-L,c) = 0;%decay_guess(r); %exchange 0 with decay guess ot make like original
                    end
                end
                shift= shift+1;
            end
            E = E_init;
            
            data_fit = zeros(length(data),1);       
            %% Fit - recovering photon arrival times and amplitudes
            startpoint = 1;
            endpoint = filtersize;       
    %         E_prime = E; not necessary
            while startpoint <= length(data)%for i = 1:length(data)        
                if endpoint > length(data)
                    endpoint = length(data);
                end  
                data_prime = data(startpoint:endpoint,1); %subset of the data
                E_prime = E(1:length(data_prime),1:length(data_prime)); %subset of the E matrix
                More = 1;
                while More == 1
                    C_prime = E_prime\data_prime;%(E_prime'*E_prime)\E_prime'*data_prime; %finding the 'concentrations' 
                    Index = find(C_prime>threshold); %keep above threshold
                    NegIndex = find(C_prime<threshold); %discard below the threshold
                    if isempty(NegIndex)
                        More = 0;
                    else
                       E_prime = E_prime(:,Index);
                    end
                end
%                 C_prime = E_prime\data_prime;%(E_prime'*E_prime)\E_prime'*data_prime; %finding the 'concentrations' 
%                 Index = find(C_prime>Threshold); %keep above threshold
%                 E_prime = E_prime(:,Index);
                conc_prime = E_prime\data_prime; %(E_prime'*E_prime)\E_prime'*data_prime;           
        %% Puts the photon amplitudes together with the arrival times
                Ampl_prime = zeros(filtersize,1);
                for j = 1:size(E_prime,2)
                    index = find(E_prime(:,j)==decay_guess(1));
                        Ampl_prime(index)= conc_prime(j);
                end
                if isempty(E_prime)
                    %disp('ITS EMPTY');
                else 
                    data(startpoint:endpoint) = data(startpoint:endpoint) - (E_prime(:,1)*Ampl_prime(1));
                end
                data_fit(startpoint,1) = Ampl_prime(1); %creates the final array of recovered photon amplitudes and arrival times
                new_start = find(Ampl_prime ~= 0);
                new_start = new_start(new_start(:) >= 2);
                if isempty(new_start) == 1 || new_start(end) <= 1
                    startpoint = startpoint + filtersize;
                    endpoint = endpoint + filtersize;
                else
                    startpoint = new_start(1) + startpoint - 1 ;
                    endpoint = endpoint + new_start(1) - 1;
                end
            end

            I = eye(filtersize,filtersize); %identity matrix for calculating the weighted convolution matrix
            M = conv2(data_fit,I);
            M = M(1:length(M)-(length(I)-1),:); %needs to be cut down to be the proper length
            %% calculates the irf and decides if the fit is good enough
            irf = M\Reduced_Data;%(M'*M)\M'*Reduced_Data;
            irf = irf/sum(irf);%(1/max(irf))*irf;
            Index2 = find(irf<0); %find negative values in the recovered irf
            irf(Index2) = 0; %set negative values to 0 
            
            if counter >= 1
                Past_Diff = Diff;
                Past_nonzero = nonzero;
                Past_X2 = X2;
            end
            Diff = norm(irf - decay_guess); %does this need to be changed to actually be the X2 value? 
            X2 = norm(Reduced_Data-data_fit);
            nonzero = length(find(data_fit>0));

            if counter >= 1
                comparison_Diff = Diff - Past_Diff;
                comparison_nonzero = nonzero - Past_nonzero;
                comparison_X2 = X2-Past_X2;
                if comparison_Diff < 0 
                    condition1 = 1;
                else
                    condition1 = 0;
                end
                if comparison_nonzero < 0 
                    condition2 = 1;
                else
                    condition2 = 0;
                end
            end
            if condition1 == 1 %&& condition2 == 1
                decay_guess = irf;
            else
%                 keepgoing = 0;
                final_Diff = Diff; %should I use this?
                final_nonzero = nonzero;
                final_X2 = Past_X2;
                [irf2,final_nonzero2,data_fit_final] = NNMF_IRF_testing_part2(irf,Reduced_Data,threshold);
                if final_nonzero2 < 1500 || sum(isnan(irf2)) > 0 
                    if sum(isnan(irf2)) > 0
                        decay_guess = IRFs(k,:)'/5 + normrnd(0,1,100,1);
                        keepgoing = 1;
                    else
                        keepgoing = 1;
                        decay_guess = irf2/5 + normrnd(0,1,100,1);
                    end
                else
                    keepgoing = 0;
                end
            end
            counter = counter + 1
        end
        results(:,k) = vertcat(MuGuess,final_Diff,final_X2,final_nonzero2);
        data_results(:,k) = data_fit_final;
        irf_results(:,k) = irf2;
        telapsed(:,k) = toc(tstart);
    end
    minimums = min(results,[],2);
    location = find(results(2,:) == minimums(2));
    recovered_events{:,n} = data_results(:,location);
    recovered_irf{:,n} = irf_results(:,location);
end

if plot_irf == 1
    for i = 1:35
    figure(i)
    plot(recovered_irf{i})
    end
end
if plot_events == 1
    for i = 1:Buffers
        if isempty(recovered_events{i})
        else
            figure(i)
            plot(temp2(value1(i)-1000:value2(i),i))
            hold on
            bar(-recovered_events{i})
            hold off
        end
    end
end

toc
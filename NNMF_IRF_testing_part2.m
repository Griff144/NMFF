%This program is to be used with the NMFF program to get the final 
function [output,output2,output3] = NNMF_IRF_testing_part2(decay_guess,Reduced_Data,Threshold)

% PerformAutocorrelation = 0;
% Reducing = 0;
% datatype = 1; %0 for .bin and 1 for xls
% plot_irf = 0;
% plot_events = 0;

% Threshold = 0.45;%4.8;%0.8;%0.0055;

%% Load data
% if datatype == 1
%     Filename = '\\10.164.16.153\Data\Scott  Griffin\deconvolution\ROC plot\TimeTraces&GroundTruth\TimeTrace2.csv';
%     Rescaled_Data = csvread(Filename);
%     Filename2 = '\\10.164.16.153\Data\Scott  Griffin\deconvolution\ROC plot\true_irf.csv';
%     GroundTruthIRF = csvread(Filename2,0,0);
%     %Filename3 = 'C:\Users\scott\OneDrive\Documents\triboluminescence\deconvolution\PeakShape.csv';
%     %PeakShape = csvread(Filename3,0,0);
    filter = 50;%100;
    reduction = 1;
%     Buffers = 1;
%     temp2 = Rescaled_Data;
% end
% irf_stru ct = load('\\10.164.16.153\Data\CaseySmith\Deconvolution\IRFs.mat');
% IRFs = irf_struct([1]).outputIrf;
% 
% recovered_events = cell(1,Buffers);
% recovered_irf = cell(1,Buffers);
% for n= 1:Buffers
% %     if value1(n) == 0
% %         n = n + 1;
% %     end
%     Reduced_Data = temp2;%temp2(value1(n)-1000:value2(n)+1000,n);
% %% Create baseline correcting high-pass filter and apply it
%     x_hp = linspace(-3,3);
%     pdf = normpdf(x_hp,0,1);
%     pdf_normalized = -pdf/sum(pdf);
%     impulse = zeros(1,100);
%     for i = 1:length(x_hp)
%         if i == 50
%             impulse(1,i) = abs(sum(pdf_normalized));
%         else
%             impulse(1,i) = 0;
%         end
%     end
%     final_filter = pdf_normalized + impulse;
%     Reduced_Data = conv(final_filter,Reduced_Data);
%     Reduced_Data = Reduced_Data(50:length(Reduced_Data)-50);
%     %% Fit 
    filtersize = filter/reduction; %size of data transient 
%     irf_results = zeros(filtersize,1);
    delta = 0;
    how_many = 1;
%     x = linspace(.1,5,filtersize);
%     %data_results = zeros(length(Reduced_Data),how_many);
    MuGuess = 1; %where the guess starts
%     % decay_guess = [0.471917366689490;0.444352660347743;0.400492898678332;0.373847308960430;0.329836788211329;0.320300677845608;0.282115586036927;0.260809998101036;0.250476052991309;0.268441947775254;0.215189596833900;0.183447499700932;0.183234410057353;0.151668624681817;0.138561988496781;0.136344574734539;0.115819953813804;0.100576173821636;0.0858093874365866;0.0976735883040951;0.0634813791347510;0.0641695240580935;0.0665161037050845;0.0468127887914506;0.0494991575730799;0.0407929310644919;0.0340367268248917;0.0374619025173877;0.0232705217253300;0.0417237091368911];
%     % decay_guess = [0.473797267322844;0.441566064852616;0.399091628639536;0.376779944700209;0.326180115998706;0.322155379691503;0.281952718949510;0.261826837967303;0.246691071297903;0.269164880138190;0.218150087392597;0.180434881277944;0.184455695811541;0.153703641120796;0.136112174938992;0.138957821489894;0.114424320327840;0.103126751622742;0.0821120016149100;0.101024099650340;0.0632004576080732;0.0614089893830797;0.0688111851444227;0.0478470450672085;0.0471602690291982;0.0439049538346425;0.0332958943075685;0.0417165023338694;0.0158823144206397;0.0515611761228003]';
%     % decay_guess = [0.473797267322844,0.441566064852616,0.399091628639536,0.376779944700209,0.326180115998706,0.322155379691503,0.281952718949510,0.261826837967303,0.246691071297903,0.269164880138190,0.218150087392597,0.180434881277944,0.184455695811541,0.153703641120796,0.136112174938992,0.138957821489894,0.114424320327840,0.103126751622742,0.0821120016149100,0.101024099650340,0.0632004576080732,0.0614089893830797,0.0688111851444227,0.0478470450672085,0.0471602690291982,0.0439049538346425,0.0332958943075685,0.0417165023338694,0.0158823144206397,0.0515611761228003]';

    for k = 1:how_many
        MuGuess = MuGuess + delta; %original guess value for the eponential pdf
        comb = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]';%1,0,0,0,0,0,0,0,0,0]';
    %    decay_guess = IRFs(k,:)';%comb .* exppdf(x,MuGuess)'; %original guess pdf 
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
                    Index = find(C_prime>Threshold); %keep above threshold
                    NegIndex = find(C_prime<Threshold); %discard below the threshold
                    if isempty(NegIndex)
                        More = 0;
                    else
                       E_prime = E_prime(:,Index);
                    end
                end
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
                keepgoing = 0;
                final_Diff = Past_Diff;
                final_nonzero = nonzero;
                final_X2 = Past_X2;
            end
            counter = counter + 1
        end
    end
% end
output = irf;
output2 = final_nonzero;
output3 = data_fit;

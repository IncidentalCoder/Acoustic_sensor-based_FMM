% Program for detecting probe movement from accelerometer data
%
clear; clc;
%
% ~~~~~~~~~~~~~~~~~~~~ LOADING THE INPUT DATA FILE ~~~~~~~~~~~~~~~~~~~~~~~%
%
% Locate and load the data files
[data_file_names,path_name] = uigetfile('*.mat','Select the data files','MultiSelect','on'); % Returns file names with the extension in a cell, and path name
addpath(path_name); % Adds the path to the file
%
if ischar(data_file_names)
    n_data_files = 1;
else
    n_data_files = length(data_file_names);
end
%%

% ~~~~~~~~~~~~~~~ INITIALIZAITON OF DETECTION VARIABLES ~~~~~~~~~~~~~~~~~~%
% VariableS for SAMPLING FREQUENCY
SF_Acstc       = 1024; % Sampling frequency of the acoustic data in Hz
SF_IMU         = 64; % Sampling frequency of the IMU data in Hz
SF_M_sensation = 1; % Sampling frequency of the Maternal sensation data in Hz
%

% Variable for thresholding
%thrsld_value_Acstc = 0.025; % Selected based on the trial-and-error check
thrsld_value_S_IMU = 0.015; % Selected based on the trial-and-error check
%

% Variables for dilation
dilation_period_S_IMU = 3.0; % dilation period in second
dilation_size_S_IMU   = round(dilation_period_S_IMU*SF_IMU); % dilation size in sample number
dilation_period_Acstc = 3.0; % dilation period in second
dilation_size_Acstc   = round(dilation_period_Acstc*SF_Acstc);
%

% Variable for matching with maternal sensation
ext_backward = 5; % Backward extension length in second for matching with maternal sensation
ext_forward = 2; % Forward extension length in second for matching with maternal sensation
reqrd_overlap = 0.3; % Matching will be considered successfull detection if the overlap (%) > this value x 100 (%)
%

% Variables for detection stats
TPD_indv = zeros(n_data_files,1);
FPD_indv = zeros(n_data_files,1);
TND_indv = zeros(n_data_files,1);
FND_indv = zeros(n_data_files,1);

DF_length = zeros(n_data_files,1); % Variable to hold data file lengths in seconds

% Loop across threshold values
ROC_analysis  = 0; % if 0, ROC analysis is not performed. Otehrwise ROC analysis is performed
thd_increment = 0.001;

if ROC_analysis
    n_loop = 100;
    thrsld_value_Acstc = 0.001;
else
    n_loop = 1;
    thrsld_value_Acstc = 0.025;
end

TPD_total = zeros(n_loop,1);
FPD_total = zeros(n_loop,1);
TND_total = zeros(n_loop,1);
FND_total = zeros(n_loop,1);

for loop_var = 1:n_loop

    fprintf('\n\tIteration: %.0f/%.0f',loop_var,n_loop) % Displaying the current value of the multiplier on screen     
    
    % Loop across the data files
    for i = 1:n_data_files

        % Loading of data files
        if n_data_files == 1
            load(data_file_names); % A cell matrix named data_var will be loaded, which contains all the data
        else
            load(data_file_names{i}); % A cell matrix named data_var will be loaded, which contains all the data
        end
        %

        % Assigning of the variables
        Acstc_data       = Subject_data.Acoustic;
        IMU_data         = Subject_data.IMU;
        M_sensation_data = Subject_data.Event;
        %
        DF_length(i) = length(Acstc_data)/SF_Acstc;

        % Removal of last 5 sec of data
        Removal_period   = 5; % in second
        Acstc_data       = Acstc_data(1:end-Removal_period*SF_Acstc,:);
        IMU_data         = IMU_data(1:end-Removal_period*SF_IMU,:);
        M_sensation_data = M_sensation_data(1:end-Removal_period*SF_M_sensation,:);
        %

        % Calculating the eucledian norm of the IMU data
        S_IMU_data = sqrt(sum(IMU_data.^2,2));
        %

        % Creation of the Maternal sensation event data variable
        M_event = find(M_sensation_data); % Finding the indices (which is the event timing in second) of the non-zero elements in M_sensation
        M_event = [M_event,ones(length(M_event),1)]; % Adding a column with a constant value of 1
        %
        [length_Acstc_data, N_Acstc_channel]= size(Acstc_data); % length and width of Acoustic data
        %

        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FILTERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %
        lowCutoff = 1;
        highCutoff = 30;
        filter_order = 4;
        %
        [b,a] = butter(filter_order,[lowCutoff highCutoff]/(SF_Acstc/2),'bandpass');
        Acstc_data_fltd= filtfilt(b,a,Acstc_data); % Zero-phase filtering
        %
        [b,a] = butter(filter_order,[1 10]/(SF_IMU/2),'bandpass');
        S_IMU_data_fltd= filtfilt(b,a,S_IMU_data); % Zero-phase filtering
        %

        % ~~~~~~~~~~~~~~~~~ EXCLUSION OF MATERNAL MOVEMENT ~~~~~~~~~~~~~~~~~~~~~~~%
        %
        % Thresholding the IMU data
        S_IMU_data_fltd_abs = abs(S_IMU_data_fltd);
        S_IMU_data_thrsld = (S_IMU_data_fltd_abs >= thrsld_value_S_IMU); % Logical column vector
        %

        % Dilating the IMU data
        SE_S_IMU = strel('line',dilation_size_S_IMU,90); % Creates a linear structuring element (vertical, as deg=90) that will have values 1 with a length of dilation_size;
        S_IMU_thd_dilated = imdilate(S_IMU_data_thrsld,SE_S_IMU); % Dilate or expand the ROI's(points with value=1) by dilation_size (half above and half below),as defined by SE
        exclusion_mask_S_IMU = imresize(double(S_IMU_thd_dilated),[length(Acstc_data_fltd),1])>0.5;
        % IMU_exclusion_mask = logical column vector with same length as the fltd_Acstc_data
        % imresize(A,[numrows numcols]) returns image that has the number of rows and columns specified by the two-element vector
        % By default, imresize uses bicubic interpolation. This is to make the mask size same as the acoustic signal size
        %
        % Application of the musk on the acoustic signal
        Acstc_data_fltd_IMU_excld = Acstc_data_fltd.*(1-exclusion_mask_S_IMU); % Maternal movement excluded data
        %

        % ~~~~~~~~~~~~~~~~~ SEGMENTATION OF THE ACOUSTIC DATA ~~~~~~~~~~~~~~~~~~~~%
        %
        % Thresholding the acoustic signal
        Acstc_data_fltd_S_IMU_excld_abs = abs(Acstc_data_fltd_IMU_excld);
        Acstc_data_thrsld = (Acstc_data_fltd_S_IMU_excld_abs >= thrsld_value_Acstc); % thresholding the acoustic data
        %
        % Combining all the thresholded acoustic colummns
        Acstc_data_thrsld_cmbd = sum(Acstc_data_thrsld,2) > 0; % Combining all the acoustic colummn
        %
        % Dilation of the thresholded acoustic data to join the near-by non-zero value
        SE_Acstc = strel('line',dilation_size_Acstc,90);
        Acstc_data_thrsld_cmbd_dilated = imdilate(Acstc_data_thrsld_cmbd,SE_Acstc);
        %
        % Labeling of dilated signal
        Acstc_data_thrsld_cmbd_dilated_labeled = bwlabel(Acstc_data_thrsld_cmbd_dilated); % Connected non-zero values are labeled as 1,2,3 etc. sequentially
        %

        % ~~~~~~~~~~~~~~~~~ MATCHING WITH MATERNAL SENSATION ~~~~~~~~~~~~~~~~~~~~~%
        %
        % Creation of the Maternal sensation map
        M_sensation_Map = zeros(size(Acstc_data_thrsld_cmbd_dilated)); % Initialization
        N_True_positive_detection = 0; % Initialization of the variable for true positive detection by the sensors
        %
        for j = 1:length(M_event) % M_event contains all the maternal sensation timings
            L   = M_event(j)*SF_Acstc; % Sample no. corresponding to the maternal sensation
            DLB = round(ext_backward*SF_Acstc); % backward extension length
            DLF = round(ext_forward*SF_Acstc); % forward extension length
            L1  = L-DLB; % sample no. for the starting point of this sensation in the map
            L2  = L+DLF; % sample no. for the ending point of this sensation in the map
            L1  = max(L1,1); % Just a check so that L1 remains higher than 1st data sample
            L2  = min(L2,length(M_sensation_Map)); % Just a check so that L2 remains lower than last data sample
            %

            % Single vector with all the sensation data mapping
            M_sensation_Map(L1:L2) = 1; % Assigns 1 to the location defined by L1:L2
            %

            % Individual checking of each maternal detection with the segmented data
            indv_M_sensation_Map        = zeros(size(Acstc_data_thrsld_cmbd_dilated)); % Initialization in every loop to remove the previous data
            indv_M_sensation_Map(L1:L2) = 1; % Only the current Maternal detection is non-zero
            curnt_matched_vector        = Acstc_data_thrsld_cmbd_dilated_labeled.*indv_M_sensation_Map; % Non-zero elements gives the matching
            curnt_matched_label         = unique(curnt_matched_vector); % Gives the label of the matched acoustic data segment

            if length(curnt_matched_label)>1 % This means, matching has been found
                curnt_matched_label = curnt_matched_label(2:end); % First element is 0, which is excluded

                for k = 1:length(curnt_matched_label) % Loop for each matched segment
                    matched_signal_original_length = length(find(Acstc_data_thrsld_cmbd_dilated_labeled == curnt_matched_label(k))); % Length of the whole matched segment
                    matched_signal_matched_length  = length(find(curnt_matched_vector == curnt_matched_label(k))); % Length of the matched portion of the segment

                    condition1 = matched_signal_matched_length > matched_signal_original_length*reqrd_overlap;% more than 30% of the segment length
                    condition2 = matched_signal_matched_length > (ext_backward+ext_forward)*SF_Acstc*reqrd_overlap; % more than 30% of the window

                    if (condition1||condition2) % Condition for overlapping
                        N_True_positive_detection = N_True_positive_detection + 1; % Increment of the true detection number
                        break % If one of the segment fulfills the condition, that Maternal detection is already detected and no more segment is checked
                    end
                end
            end
        end
        %

        % ~~~~~~~~~~~~ EXTRACTION OF THE MATCHED DATA FOR T-F ANALYSYS ~~~~~~~~~~~%
        %
        % Determination of the matched vector and the total no. of matched segment
        matched_vector                   = Acstc_data_thrsld_cmbd_dilated_labeled.*M_sensation_Map; % Contains all the matched acoustic segments irrespective of the overlapping condition
        matched_original_segment_numbers = unique(matched_vector); % Index of the matched segments irrespective of the overlapping condition
        matched_original_segment_numbers = matched_original_segment_numbers(2:end); % removal of first index, which is zero
        N_total_matching                 = length(matched_original_segment_numbers);
        % This gives the matched segment number in the thresolded acoustic data
        %

        % Extraction of matched acoustic data
        extracted_all_matched_Acstc_data = Acstc_data_fltd_IMU_excld.*(matched_vector>0); % Contains all the matched data irrespective of % of overlapping
        % Holds the undilated original IMU_exclded Acstc data that matched with maternal sensation
        %
        % Checking of the percentage of matching and determining the best channel for each matched segments and putting them into seperate columns
        extraction_size_indv = 10; % defines the size (in second) of each column in the extracted data
        pre_extnd_length     = 1; % extension length for data extraction before the segment starts
        post_extnd_length    = 2; % extension length for data extraction after the segment ends
        clmn_index           = 1; % Initialization; this variable will decide in which column the extracted segment will be put on
        %
        extracted_full_segment_Acstc_data_in_seprt_colmn       = zeros(extraction_size_indv*SF_Acstc,1); % Initialization
        extracted_full_segment_Acstc_data_in_seprt_colmn_extnd = zeros(extraction_size_indv*SF_Acstc,1); % Initialization
        properly_matched_full_segment_vector                   = zeros(size(matched_vector)); % Initialization
        %
        for j = 1:N_total_matching % Looping through each matched data set
            % Indexes
            index_current_indv_matched = find(matched_vector == matched_original_segment_numbers(j)); % This has the index of the matched part of the segment only
            index_curnt_indv_matched_full_segment = find(Acstc_data_thrsld_cmbd_dilated_labeled == matched_original_segment_numbers(j));
            % This has the index of the full segment for the current matched segment
            %
            % Extension to get some part of the data before and after the segment, which makes the Spectrogram look better
            pre_extnd_index  = (max((index_curnt_indv_matched_full_segment(1)- pre_extnd_length*SF_Acstc),1):1:(index_curnt_indv_matched_full_segment(1)-1))';
            post_extnd_index = ((index_curnt_indv_matched_full_segment(end)+ 1):1:min((index_curnt_indv_matched_full_segment(end)+post_extnd_length*SF_Acstc),length_Acstc_data))';
            index_curnt_indv_matched_full_segment_extnd = [pre_extnd_index;index_curnt_indv_matched_full_segment;post_extnd_index];
            %
            % Lengths
            length_current_indv_matched                  = length(index_current_indv_matched);
            length_curnt_indv_matched_full_segment       = length(index_curnt_indv_matched_full_segment);
            length_curnt_indv_matched_full_segment_extnd = length(index_curnt_indv_matched_full_segment_extnd);
            %
            energy_previous_clmn_indv_matched_full_segment = 0; % Initializing the energy of the current set of matched data
            %
            condition1 = length_current_indv_matched > length_curnt_indv_matched_full_segment*reqrd_overlap;% more than 30% of the segment length
            condition2 = length_current_indv_matched > (ext_backward+ext_forward)*SF_Acstc*reqrd_overlap; % more than 30% of the window
            if (condition1||condition2) % if true, the condition for overlap is met
                for k = 1 : N_Acstc_channel % Looping through each column for selecting the best
                    curnt_clmn_indv_matched_full_segment_Acstc_data = Acstc_data_fltd_IMU_excld(index_curnt_indv_matched_full_segment,k);
                    % data for current full segment; energy will be calculated for full segment
                    curnt_clmn_indv_matched_full_segment_Acstc_data_extended = Acstc_data_fltd(index_curnt_indv_matched_full_segment_extnd,k);
                    % data for current full segment with extension; also IMU exclusion is not considered, just to have a good view in the spectrogram
                    properly_matched_full_segment_vector(index_curnt_indv_matched_full_segment) = 1;
                    % This vector contains the properly matched full segment data as non-zeor value (=1)
                    % This will be necessary to extract the false positive detections
                    %
                    energy_curnt_clmn_indv_matched_full_segment = sum(curnt_clmn_indv_matched_full_segment_Acstc_data.^2); % updating the energy variable
                    %
                    % Finding the best channel based on maximum energy
                    if energy_curnt_clmn_indv_matched_full_segment > energy_previous_clmn_indv_matched_full_segment % Condition for finding the best channel
                        extracted_full_segment_Acstc_data_in_seprt_colmn(1:length_curnt_indv_matched_full_segment,clmn_index) = curnt_clmn_indv_matched_full_segment_Acstc_data;
                        % Putting individual matched data in each column
                        extracted_full_segment_Acstc_data_in_seprt_colmn_extnd(1:length_curnt_indv_matched_full_segment_extnd,clmn_index) = curnt_clmn_indv_matched_full_segment_Acstc_data_extended;
                        % Putting individual matched data in each column
                        energy_previous_clmn_indv_matched_full_segment = energy_curnt_clmn_indv_matched_full_segment;
                        % Updating the energy variable
                    end
                end
                clmn_index = clmn_index + 1; % Increment of the column index for storing the next segment in the next column
            end
        end
        %
        % s = spectrogram(indv_extracted_full_segment_Acstc_data(:,1));
        % spectrogram(indv_extracted_full_segment_Acstc_data(1:length(indv_matched_full_segment_Acstc_data),1),'yaxis')
        %

        %~~~~~~~~~~~ DETERMINATION OF TRUE NEGATIVE AND FALSE POSITIVE~~~~~~~~~~~~~
        %
        index_total_negative_segments = find(M_sensation_Map == 0);
        total_thrsld_negative_segments = Acstc_data_thrsld_cmbd_dilated(index_total_negative_segments);
        % In this vector, only the regions of (true negative + false positive) detections are non-zero (= 1)
        length_total_thrsld_negative_segments = length(total_thrsld_negative_segments);
        %
        window_size = ext_backward + ext_forward;
        % window size is equal to the window size used to create the maternal sensation map
        index_window_start = 1;
        index_window_end = index_window_start + SF_Acstc*window_size;
        N_False_positive_detection = 0; % Initialization
        N_True_negative_detection = 0; % Initialization
        minm_overlap = 1.25; % Minimum overlap in second
        %
        while (index_window_start < length_total_thrsld_negative_segments)
            indv_negative_segment_window = total_thrsld_negative_segments (index_window_start: index_window_end);
            index_non_zero = find(indv_negative_segment_window>0);
            if length(index_non_zero) >= (minm_overlap*SF_Acstc)
                N_False_positive_detection = N_False_positive_detection + 1;
            else
                N_True_negative_detection = N_True_negative_detection + 1;
            end
            index_window_start = index_window_end + 1;
            index_window_end = min(index_window_start + SF_Acstc*window_size,length_total_thrsld_negative_segments);
        end
        %

        % ~~~~~~~~~~~~~~~ EXTRACTION OF FALSE POSITIVE SIGNALS ~~~~~~~~~~~~~~~~~~~%
        %
        False_positive_thrsld_vector = Acstc_data_thrsld_cmbd_dilated.*(1-properly_matched_full_segment_vector); % Binary vector with false positives
        False_positive_thd_vector_labeled = bwlabel(False_positive_thrsld_vector);
        %
        extracted_False_positive_Acstc_data = Acstc_data_fltd_IMU_excld.*False_positive_thrsld_vector; % Contains all the false positives
        %
        False_positives_labels = unique(False_positive_thd_vector_labeled);
        False_positives_labels = False_positives_labels(2:end); % Removing the first label, which is zero
        N_False_positive = length(False_positives_labels);
        %
        pre_extnd_length = 1; % extension length for data extraction before the segment starts
        post_extnd_length = 2; % extension length for data extraction after the segment ends
        %
        extracted_False_positive_Acstc_data_in_seprt_clmn = zeros(10*SF_Acstc,1); % Initialization
        extracted_False_positive_Acstc_data_in_seprt_clmn_extnd = zeros(10*SF_Acstc,1); % Initialization
        %
        for j = 1 : N_False_positive
            index_currrent_indv_False_positive = find(False_positive_thd_vector_labeled == j);
            current_indv_False_positive_Acstc_data = Acstc_data_fltd_IMU_excld(index_currrent_indv_False_positive,:);
            %
            % Extension to get some part of the data before and after the segment, which makes the Spectrogram look better
            pre_extnd_index = (max((index_currrent_indv_False_positive(1)- pre_extnd_length*SF_Acstc),1):1:(index_currrent_indv_False_positive(1)-1))';
            post_extnd_index = ((index_currrent_indv_False_positive(end)+ 1):1:min((index_currrent_indv_False_positive(end)+post_extnd_length*SF_Acstc),length(Acstc_data_fltd)))';
            index_currrent_indv_False_positive_extnd = [pre_extnd_index;index_currrent_indv_False_positive;post_extnd_index];
            current_indv_False_positive_Acstc_data_extnd = Acstc_data_fltd(index_currrent_indv_False_positive_extnd,:);
            % data for current full segment with extension; also IMU exclusion is not considered, just to have a good view in the spectrogram
            %
            % Lenght determination
            length_current_indv_False_positive_Acstc_data = length(current_indv_False_positive_Acstc_data);
            length_current_indv_False_positive_Acstc_data_extnd = length (current_indv_False_positive_Acstc_data_extnd);
            %
            energy_previous_clmn_False_positive = 0;
            %
            % Loop for finding the best channel based on maximum energy
            for k = 1:N_Acstc_channel
                current_clmn_indv_False_positive_Acstc_data = current_indv_False_positive_Acstc_data(:,k);
                current_clmn_indv_False_positive_Acstc_data_extnd = current_indv_False_positive_Acstc_data_extnd (:,k);
                %
                energy_current_clmn_False_positive = sum(current_clmn_indv_False_positive_Acstc_data.^2);
                %
                % Loop to find the channel with maximum energy
                if energy_current_clmn_False_positive > energy_previous_clmn_False_positive % Condition for finding the best channel
                    extracted_False_positive_Acstc_data_in_seprt_clmn(1:length_current_indv_False_positive_Acstc_data,j) = current_clmn_indv_False_positive_Acstc_data;
                    % Putting individual matched data in each column
                    extracted_False_positive_Acstc_data_in_seprt_clmn_extnd(1:length_current_indv_False_positive_Acstc_data_extnd,j) = current_clmn_indv_False_positive_Acstc_data_extnd;
                    % Putting individual matched data in each column
                    energy_previous_clmn_False_positive = energy_current_clmn_False_positive;
                    % Updating the energy variable
                end
            end
        end
        %

        %~~~~~~~~~~~~~~~~~~ DETERMINATION OF FALSE NEGATIVES ~~~~~~~~~~~~~~~~~%
        %
        N_total_Sensor_detection = max(bwlabel(Acstc_data_thrsld_cmbd_dilated)); % No. of detection by the sensor;
        N_Maternal_detection = length(M_event); % No of detection by maternal sensation
        %
        N_False_negative_detection = N_Maternal_detection - N_True_positive_detection;
        %

        % ~~~~~~~~~~~~~~ DETECTION STATS FOR INDIVIDUAL DATA SETS ~~~~~~~~~~~~%
        TPD_indv(i) = N_True_positive_detection;
        FPD_indv(i) = N_False_positive_detection;
        TND_indv(i) = N_True_negative_detection;
        FND_indv(i) = N_False_negative_detection;

        % Storing the extracted TPD and FPD data
        if ROC_analysis == 0            
            TPD_extracted{i} = extracted_full_segment_Acstc_data_in_seprt_colmn;
            FPD_extracted{i} = extracted_False_positive_Acstc_data_in_seprt_clmn;
        end
        %

    end
    %

    % Overal detection stats
    TPD_total(loop_var) = sum(TPD_indv);
    FPD_total(loop_var) = sum(FPD_indv);
    TND_total(loop_var) = sum(TND_indv);
    FND_total(loop_var) = sum(FND_indv);

    % Performance evaluation
    if ROC_analysis == 0       
        SEN = TPD_total/(TPD_total+FND_total);
        PPV = TPD_total/(TPD_total+FPD_total); % same as precision
        F1_score = 2*SEN*PPV/(SEN+PPV);
        SPE = TND_total/(TND_total+FPD_total);
        ACC = (TPD_total+TND_total)/(TPD_total+TND_total+FPD_total+FND_total);
        FPR = 1 - SPE; % Flase positive rate
    end
    %

    % Update the threshold value
    thrsld_value_Acstc = thrsld_value_Acstc + thd_increment;
end

fprintf('\nAnalysis completed.\n') % Displaying the current value of the multiplier on screen

%% ~~~~~~~~~~~~~~~~~~~~~ TIME-FREQUENCY ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~%
% PSD estimation: 
% Use non-extended extracted signals for PSD estimation

% Parameters and variables for PSD estimation
s_size = SF_Acstc/2;
w = hann(s_size); % Hann window of length 1s
n_overlap = s_size/2; % 50% overlap between the segments
p_overlap = 0.5; % 50% overlap between the segments
nfft = 2*s_size; % number of fourier points will be same as sampling freq. Therefore no zero padding

Pxx_TPD_avg = zeros(SF_Acstc/2 + 1, 1);
Pxx_detrended_TPD_avg = zeros(SF_Acstc/2 + 1, 1);
Pxx_FPD_avg = zeros(SF_Acstc/2 + 1, 1);
Pxx_detrended_FPD_avg = zeros(SF_Acstc/2 + 1, 1);
n_TPD_cycle = 0;
n_FPD_cycle = 0;

for i = 1:length(TPD_extracted)
    for j = 1 : size(TPD_extracted{i},2)

        data = TPD_extracted{i}(:,j);

        [Pxx_TPD, f_TPD] = pwelch(data, w, n_overlap, nfft, SF_Acstc);
        [Pxx_detrended_TPD, f_detrended_TPD] = pwelch_new(data, w, p_overlap, nfft, SF_Acstc, 'half', 'plot', 'mean');
        % pwelch_new() is a user defined function that allows detrending.
        %   Input conditions: half - one sided PSD; plot- normal plotting;
        %                     mean-  remove the mean value of each segment from each segment of the data.

        Pxx_TPD_avg = Pxx_TPD_avg + Pxx_TPD; % Summing the PSD for each TPD
        Pxx_detrended_TPD_avg = Pxx_detrended_TPD_avg + Pxx_detrended_TPD;

        n_TPD_cycle = n_TPD_cycle + 1;
    end
end

for i = 1:length(FPD_extracted)
    for j = 1 : size(FPD_extracted{i},2)

        data = FPD_extracted{i}(:,j);

        [Pxx_FPD, f_FPD] = pwelch(data, w, n_overlap, nfft, SF_Acstc);
        [Pxx_detrended_FPD, f_detrended_FPD] = pwelch_new(data, w, p_overlap, nfft, SF_Acstc, 'half', 'plot', 'mean');
        % pwelch_new() is a user defined function that allows detrending.
        %   Input conditions: half - one sided PSD; plot- normal plotting;
        %                     mean-  remove the mean value of each segment from each segment of the data.

        Pxx_FPD_avg = Pxx_FPD_avg + Pxx_FPD; % Summing the PSD for each TPD
        Pxx_detrended_FPD_avg = Pxx_detrended_FPD_avg + Pxx_detrended_FPD;

        n_FPD_cycle = n_FPD_cycle + 1;
    end
end

% Averaging the PSD summation
n_TPD_cycle           = n_TPD_cycle - 1; % One of the extracted data set did not have any signal 5th element in the cell matrix
Pxx_TPD_avg           = Pxx_TPD_avg/n_TPD_cycle; 
Pxx_detrended_TPD_avg = Pxx_detrended_TPD_avg/n_TPD_cycle;
Pxx_FPD_avg           = Pxx_FPD_avg/n_FPD_cycle;
Pxx_detrended_FPD_avg = Pxx_detrended_FPD_avg/n_FPD_cycle;

% Plotting the PSD
x_lim = 50;
subplot(2,1,1)
for i = 1:1    
    plot(f_detrended_TPD, Pxx_detrended_TPD_avg(:,i))
    xlim([0 x_lim])
    hold on;
    plot(f_TPD, Pxx_TPD_avg(:,i))
end
hold off;

subplot(2,1,2)
for i = 1:1    
    plot(f_detrended_FPD, Pxx_detrended_FPD_avg(:,i))
    xlim([0 x_lim])
    hold on;
    plot(f_FPD, Pxx_FPD_avg(:,i))
end
hold off;
%

%% Spectrogram
% Use extended extracted signals for PSD estimation

% Parameters and variables for PSD estimation
s_size = SF_Acstc/2; % Size of each segment in the STFT
w = hann(s_size); % Hann window of length 0.5 S
n_overlap = floor(s_size/1.25); % 80% overlap between the segments
nfft = s_size*2; % number of fourier points will be twice sampling freq. This pads 0 at the end to increase the frequency resolution
% For the above setting, the obtained time resolution = 100 ms, and frequency resolution = 1 Hz 

% Selection of data file and TPD
fst_data_file_no = 6;
fst_TPD_no       = 11;
snd_data_file_no = 1;
snd_TPD_no       = 1;

t_start = 0; % Starting time in s
t_end = 7; % Ending time in s

% Good candidates: 
%   TPD: (File_no, TPD_no): (6,11), (1,1) 
%   FPD: (File_no, TPD_no): (2,1), (6,8)


% Plot settings
legend_all={'Acoustic sensor'};
L_width = 4; % Width of the box
p_L_width = 5; % Width of the lines in the plot
F_size = 28; % Font size
F_name = 'Arial';

tiledlayout(2,2,'Padding','tight','TileSpacing','loose'); % 3x3 tile in most compact fitting

% Get the data to be plotted from the extracted data set
fst_TPD_data = TPD_extracted{fst_data_file_no}(:, fst_TPD_no); % TO get the data for the right sensor only
% TPD_data = fst_TPD_data(t_start*SF_Acstc:t_end*SF_Acstc+1,1); % uncomment if you want to plot within a range defined by t_start:t_end
snd_TPD_data = TPD_extracted{snd_data_file_no}(:, snd_TPD_no); % TO get the data for the right sensor only
% snd_TPD_data = snd_TPD_data(t_start*SF_Acstc:t_end*SF_Acstc+1,1); % uncomment if you want to plot within a range defined by t_start:t_end


% Plot the sensor data in the first row ---------------------------------
% first TPD data
time_vec = (0:length(fst_TPD_data)-1)/SF_Acstc;
nexttile
plot(time_vec, fst_TPD_data, 'LineWidth', p_L_width)

% Set the axis properties
ylabel('Amplitude (a.u.)'); % Ylabel is only applied to the leftmost axis
ax = gca; % To get the coordinate axis
ax.YAxis.Exponent = -2;
ylim([-0.1 0.1])

%xlabel('Time(s)');
xlim([0 6])
%xticks(0:1:10);
set(gca, 'FontName', F_name, 'FontSize', F_size, 'linewidth', L_width) % Sets the font type and size of the labels

% second TPD data
time_vec = (0:length(snd_TPD_data)-1)/SF_Acstc;
nexttile
plot(time_vec, snd_TPD_data, 'LineWidth', p_L_width)

% Set the axis properties
%ylabel('Amplitude (a.u.)'); % Ylabel is only applied to the leftmost axis
ax = gca; % To get the coordinate axis
ax.YAxis.Exponent = -2;
ylim([-0.04 0.04])

%xlabel('Time(s)');
xlim([0 6])
%xticks(0:1:10);
set(gca, 'FontName', F_name, 'FontSize', F_size, 'linewidth', L_width) % Sets the font type and size of the labels


% Plot the spectrogram in the seond row --------------------------------
[~,f,t,P] = spectrogram(fst_TPD_data, w, n_overlap, nfft, SF_Acstc, 'psd', 'yaxis'); % P is the PSD
nexttile
% spectrogram(data,w_spec,n_overlap_spec,nfft_spec,Fs_sensor,'yaxis') % Plot using spectrogram
imagesc(t, f, (P+eps)) % Plot using imagesc.Add eps like pspectrogram does: eps is the distance between 1.0 and the next largest double precision number (eps = 2.2204e-16)
axis xy % Corrects the Y-axis order: By default imagesc() puts Y-axis to high-to-low order
h = colorbar; % To display the colorbar on the left

% Set the axis properties
%h.Label.String = 'PSD (a.u.^2/Hz)'; % Labels the colorbar
%h.Ticks = 3*10^-3:3*10^-3:12*10^-3  ; % To manage colorbar ticks

colormap(jet)
ylim([0,30]);
%yticks(0:10:50);

xlim([-inf 6])
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'FontName', F_name, 'FontSize', F_size, 'linewidth', L_width) % Sets the font type & size, and the line width
%

% 2nd Spectrogram
[~,f,t,P] = spectrogram(snd_TPD_data, w, n_overlap, nfft, SF_Acstc, 'psd', 'yaxis'); % P is the PSD
nexttile
% spectrogram(data,w_spec,n_overlap_spec,nfft_spec,Fs_sensor,'yaxis') % Plot using spectrogram
imagesc(t, f, (P+eps)) % Plot using imagesc.Add eps like pspectrogram does: eps is the distance between 1.0 and the next largest double precision number (eps = 2.2204e-16)
axis xy % Corrects the Y-axis order: By default imagesc() puts Y-axis to high-to-low order
h = colorbar; % To display the colorbar on the left

% Set the axis properties
h.Label.String = 'PSD (a.u.^2/Hz)'; % Labels the colorbar
%h.Ticks = 3*10^-3:3*10^-3:12*10^-3  ; % To manage colorbar ticks

colormap(jet)
ylim([0,30]);
%yticks(0:10:50);

xlim([-inf 6])
xlabel('Time (s)');
%ylabel('Frequency (Hz)');
set(gca, 'FontName', F_name, 'FontSize', F_size, 'linewidth', L_width) % Sets the font type & size, and the line width
%

%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA SAVING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% %
% filename = input('Save data as: ', 's');
% % save (filename,'extracted_full_segment_Acstc_data_in_seprt_colmn', '-ascii');
% % save (filename,'extracted_full_segment_Acstc_data_in_seprt_colmn_extnd', '-ascii');
% save (filename,'extracted_False_positive_Acstc_data_in_seprt_clmn', '-ascii');
% % save (filename,'extracted_False_positive_Acstc_data_in_seprt_clmn_extnd', '-ascii');
% %
%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOTTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
% Time vector definition
time_adc = (1/SF_Acstc:1/SF_Acstc:(length(Acstc_data)/SF_Acstc))';
time_imu = (1/SF_IMU:1/SF_IMU:(length(IMU_data)/SF_IMU))';
%
% Subplot
figure(1)
subplot(4,1,1)
plot(time_imu,S_IMU_data_fltd,time_imu,S_IMU_thd_dilated*0.05, 'r');
legend('IMU data','Maternal movement');
legend off;
grid minor;
ylim([0 0.06]);
ylabel ('IMU Signal');
%
subplot(4,1,2)
p1 = plot(time_adc, Acstc_data_fltd(:,1),time_adc, Acstc_data_fltd(:,2),time_adc, Acstc_data_fltd(:,3),...
    time_adc, Acstc_data_fltd(:,4),'DisplayName','Acstc Data');
legend off;
hold on;
p2 = plot (time_adc,exclusion_mask_S_IMU*0.05,'r',...
    M_event(:,1),M_event(:,2)*0.04,'mo','MarkerSize',4 );
legend(p2,{'IMU exclusion mask','Maternal detection'});
legend off;
ylim([0 0.08]);
grid minor;
ylabel ('Acoustic signal');
hold off;
%
subplot(4,1,3)
p1 = plot(time_adc, Acstc_data_fltd_IMU_excld(:,1),time_adc, Acstc_data_fltd_IMU_excld(:,2),...
    time_adc, Acstc_data_fltd_IMU_excld(:,3),time_adc, Acstc_data_fltd_IMU_excld(:,4),...
    'DisplayName','IMU excld Acstc data');
legend off;
hold on;
p2 = plot(time_adc,Acstc_data_thrsld_cmbd_dilated*0.05,time_adc,matched_vector*0.06,time_adc,M_sensation_Map*0.07,...
    M_event(:,1),M_event(:,2)*0.04,'mo','MarkerSize',4);
legend(p2,{'Segmntd acstc','Matching vector','Senstn Map','Matrnl dtctn'});
legend off;
ylim([0 0.08]);
grid minor;
hold off;
%
subplot(4,1,4)
p1 = plot(time_adc, extracted_all_matched_Acstc_data, 'DisplayName', 'Extracted data');
legend off;
hold on;
p2 = plot(M_event(:,1),M_event(:,2)*0.04,'mo','MarkerSize',4);
legend(p2,'Maternal detection');
legend off;
ylim([0 0.08]);
grid minor;
ylabel ('Fetal Movement');
xlabel ('Time(s)')
plotbrowser('on');
%
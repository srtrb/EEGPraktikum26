%% EEG PREPROCESSING PIPELINE
%{
-----------------------------------------------------------------------------------------------
Enter your project root ('data') and the spm path in the Set path part -> XXX !

Project folder strcuture should look like this:
    data /
        00Behavioural /
            images
            logs
            neuronavigation /
                Gian_ID04.sfp
        01EEG /
            raw /
                SPNCartoons_ID04.bdf
            spm /
                channelselection.mat
                avref_eog.mat
                trialdef.mat
        05Anat

When asked to enter channels for interpolation, enter: {'P5', 'P10'}
------------------------------------------------------------------------------------------------
%}
clear; clc; 

% Set path
project_root = ['C:\XXX'];
spm_path = ['C:\XXX'];
addpath(spm_path);
spm('defaults', 'EEG');
addpath(fullfile(spm_path, 'external', 'fieldtrip'));
cd(project_root);

%% 0. Function for interpolating bad channels
% --------------------------------------------------------------------

function D2 = spm_interpolate_bad_channels(D)
% SPM_INTERPOLATE_BAD_CHANNELS - Load SPM EEG file, let user mark bad channels,
% interpolate them using spline (uses embedded sensors location in meeg file), and save the updated data.
%
% Inputs:
%   D    - meeg file loaded in workspace
% Outputs:
%   D    - meeg file with specified channels interpolated (also saved on
%   disk in same folder as loaded D file with prefix "interpolate_"

    data = spm2fieldtrip(D);

    cfg = [];
    cfg.length = 10;
    cfg.overlap = 0;
    %data_epoched = ft_redefinetrial(cfg, data);
    data_epoched = data;
    cfg = [];
    cfg.preproc.demean = 'yes';
    cfg.preproc.lpfilter = 'yes'; 
    cfg.preproc.lpfreq = 45; 
    cfg.preproc.hpfilter = 'yes'; 
    cfg.preproc.hpfreq = 1; 
    cfg.preproc.hpinstabilityfix = 'reduce'; 
    cfg.ylim = [-20 20];
    % ft_databrowser(cfg, data_epoched);
    if isfield(cfg,'colormap')
        cfg = rmfield(cfg,'colormap');
    end  % optional safety
    % ft_databrowser(cfg, data_epoched);  % COMMENT OUT

    % Let user input bad channels
    disp('Channel labels:');
    bad_labels = input('Enter bad channels as a cell array (e.g. {''F3'', ''T7''}): ');

    for i =1:length(bad_labels)
        if any(strcmp(bad_labels{i}, data.label)) == 0
            error(sprintf('The typed channel: %s do not exist', bad_labels{i}))
        end
    end

    if length(bad_labels) > 0
        cfg               = [];
        cfg.method = 'spline';
        cfg.badchannel    = bad_labels;
        %cfg.neighbours = neighbours;
        data_corr = ft_channelrepair(cfg, data);
       
    else
        data_corr = data;
    end

    D2 = D.copy(['interpolate_' fname(D)]);

    %check that dimension match
    if numel(indchantype(D, 'EEG')) == numel(data_corr.label)
        D2(indchantype(D,'EEG'),:) = data_corr.trial{1,1};
        D2.save();
    else
        error('something went wrong with channel indices')
    end
    
    fprintf('Done. Saved interpolated data as: %s\n', fullfile(D2.path, [D2.fname]));
end


%% 0. Function for classifying trials
% --------------------------------------------------------------------

function [new_trl, new_conditionlabels, condition_table] = classify_roving_trials(trl, conditionlabels)
% CLASSIFY_ROVING_TRIALS
%
% Classifies trials from a roving oddball paradigm into:
% For each rove:
%   - LAST trial = STANDARD
%   - FIRST trial of the following rove = DEVIANT
%
% Rove length is classified as:
%   n2
%   n3
%   n4-5
%   n6-8
%   n9plus
%
% For STANDARD trials:
% n = length of the current rove
%
% For DEVIANT trials:
% n = length of the preceding rove
% 
% The first trial of the experiment cannot be a deviant because there is no preceding rove.
%
% INPUTS
%   trl               Trial definition matrix returned by spm_eeg_definetrial
%   conditionlabels   condition labels cell array returned by spm_eeg_definetrial
%
% OUTPUTS
%   new_trl
%   new_conditionlabels
%   condition_table

% Convert condition labels to numeric values
values = nan(length(conditionlabels),1);

values(strcmp(conditionlabels,'low_intensity'))  = 1;
values(strcmp(conditionlabels,'high_intensity')) = 2;

if any(isnan(values))
    error('Some condition labels are neither low_intensity nor high_intensity.');
end

n = length(values);

% Find rove boundaries
% A new rove starts whenever the stimulation intensity changes 
rove_start = [1; find(diff(values) ~= 0) + 1]; 

% The end of each rove is the trial immediately before the next rove 
rove_end = [rove_start(2:end)-1; n]; 

n_roves = length(rove_start);

% Initialise output variables 
keep_idx = []; 
new_conditionlabels = {}; 

condition = {}; 
intensity = {}; 
rove_length = []; 
rove_number = [];

% Loop through roves 
for r = 1:n_roves 
    start_idx = rove_start(r); 
    end_idx = rove_end(r); 
    
    % Number of stimuli in this rove 
    current_rove_length = end_idx - start_idx + 1; 
    
    % Intensity 
    if values(start_idx) == 1 
        current_intensity = 'low'; 
    else 
        current_intensity = 'high'; 
    end

    % Convert rove length to n2 / n3 ...
    if current_rove_length == 2 
        n_label = 'n2'; 
    elseif current_rove_length == 3 
        n_label = 'n3'; 
    elseif current_rove_length == 4 
        n_label = 'n4-5'; 
    elseif current_rove_length == 5 
        n_label = 'n4-5'; 
    elseif current_rove_length == 6 
        n_label = 'n6-8'; 
    elseif current_rove_length == 7 
        n_label = 'n6-8'; 
    elseif current_rove_length == 8 
        n_label = 'n6-8'; 
    elseif current_rove_length >= 9 
        n_label = 'n9plus'; 
    
    else 
        % Rove of length 1 
        n_label = ''; 
    end

    % LAST TRIAL OF ROVE = STANDARD
    % A standard requires at least 2 repetitions in the rove 
    
    if current_rove_length >= 2 
        trial_idx = end_idx; 
    
        keep_idx(end+1) = trial_idx; 

        % SPM condition label
        new_conditionlabels{end+1,1} = sprintf( ...
            'standard_%s_%s', ...
            current_intensity, ...
            n_label);

        % Metadata
        condition{end+1,1} = 'standard'; 
        intensity{end+1,1} = current_intensity; 
        rove_length(end+1,1) = current_rove_length; 
        rove_number(end+1,1) = r; 
    end

    % FIRST TRIAL OF NEXT ROVE = DEVIANT 
    
    if r < n_roves  
    
        deviant_idx = rove_start(r+1); 
    
        if values(deviant_idx) == 1 
            deviant_intensity = 'low'; 
        else 
            deviant_intensity = 'high'; 
        end 
        
        preceding_length = current_rove_length;

        if preceding_length == 2
            deviant_n_label = 'n2';
        elseif preceding_length == 3
            deviant_n_label = 'n3';
        elseif preceding_length == 4
            deviant_n_label = 'n4-5';
        elseif preceding_length == 5
            deviant_n_label = 'n4-5';
        elseif preceding_length == 6
            deviant_n_label = 'n6-8';
        elseif preceding_length == 7
            deviant_n_label = 'n6-8';
        elseif preceding_length == 8
            deviant_n_label = 'n6-8';
        elseif preceding_length >= 9
            deviant_n_label = 'n9plus';
            
        else
            deviant_n_label = '';
        end
           
        if preceding_length >= 2

             keep_idx(end+1) = deviant_idx; 

            % SPM condition labels
            new_conditionlabels{end+1,1} = sprintf( ...
                'deviant_%s_%s', ...
                deviant_intensity, ...
                deviant_n_label);

            % Metadata
            condition{end+1,1} = 'deviant'; 
            intensity{end+1,1} = deviant_intensity;
            rove_length(end+1,1) = preceding_length;
            rove_number(end+1,1) = r+1;
        
        end

    end

end

% Sort retained trials chronologically 
[keep_idx, sort_order] = sort(keep_idx); 

new_conditionlabels = new_conditionlabels(sort_order); 

condition = condition(sort_order); 
intensity = intensity(sort_order); 
rove_length = rove_length(sort_order); 
rove_number = rove_number(sort_order);

% Create trial definitions 
new_trl = trl(keep_idx,:); 

% Create clean condition table 
condition_table = table( ... 
    keep_idx(:), ... 
    new_conditionlabels(:), ...
    condition(:), ...
    intensity(:), ...
    rove_length(:), ... 
    rove_number(:), ... 
    'VariableNames', { ... 
        'trial_index', ... 
        'spm_condition', ...
        'condition', ... 
        'intensity', ... 
        'rove_length', ... 
        'rove_number'}); 

% Display classification 

fprintf('\nRoving oddball classification\n'); 
fprintf('====================================================\n'); 

disp(condition_table); 

fprintf('====================================================\n'); 
fprintf('Original trials : %d\n', n); 
fprintf('Retained trials : %d\n', length(keep_idx)); 
fprintf('Excluded trials : %d\n\n', n-length(keep_idx)); 
fprintf('\nSPM conditions:\n');

unique_conditions = unique(new_conditionlabels, 'stable');

for i = 1:length(unique_conditions)
    fprintf('   %s\n', unique_conditions{i});
end

fprintf('\n');
end



%% 1. Continuous preprocessing
% --------------------------------------------------------------------

%convert
S = []; 
S.dataset = fullfile(project_root, '01EEG', 'raw', 'SPNCartoons_ID04.bdf');
D = spm_eeg_convert(S);

%display_SPM_data(D)

% select channels (added now for completeness)
load(fullfile(project_root, 'channelselection.mat'));
S = []; 
S.D = D; 
S.channels = label;
S.prefix = 'c'; 
D = spm_eeg_crop(S); 

D = chantype(D, 65, 'EOG');
D = chantype(D, 66, 'EOG');
D = chantype(D, 67, 'EOG');
D = chantype(D, 68, 'EOG');
D.save(); 

%display_SPM_data(D)

% Prepare (load sensor file)
S = [];
S.D = D;

S.task = 'loadeegsens';
S.source = 'locfile';
S.sensfile = fullfile(project_root, '00Behavioural', 'neuronavigation', 'Gian_ID04.sfp');
S.prefix = 'p'; 
D = spm_eeg_prep(S);


%filter (HP)
S = []; 
S.D = D; 
S.band = 'high'; 
S.freq = 0.01; 
D = spm_eeg_filter(S); 

%display_SPM_data(D)


%downsample
S = [];
S.D = D;
S.fsample_new = 200; 
D = spm_eeg_downsample(S); 

%display_SPM_data(D)

%% NEW PART - INTERPOLATE BAD CHANNELS

D = spm_interpolate_bad_channels(D);

%display_SPM_data(D)

%% END - INTERPOLATE BAD CHANNELS


%% Prepare, montage / re-referencing

load(fullfile(project_root, 'avref_eog.mat'));
S = []; 
S.D = D; 
S.montage = montage; 
D = spm_eeg_montage(S); 

%display_SPM_data(D)


%filter(LP)
S = []; 
S.D = D; 
S.band = 'low'; 
S.freq = 48; 
D = spm_eeg_filter(S); 

%display_SPM_data(D)



%% 2. Eye blink correction on continuous data 
% --------------------------------------------------------------------

% Ensure SPM is in the path
S = [];
S.D = D;
S.mode = 'mark'; % Change 'Mode' to 'Mark'
S.methods.fun = 'eyeblink'; % Detection algorithm
S.methods.settings.threshold = 4;
S.methods.channels = 'HEOG';
S.methods.settings.excwin = 0;
D_ebf = spm_eeg_artefact(S);

%check the events that have been added to the file
display_SPM_data(D_ebf)

% create epoched events around eyeblinks (and then average them)
S = []; 
S.D = D_ebf; 
S.timewin = [-500 500];
S.trialdef(1).conditionlabel = 'Eyeblink'; 
S.trialdef(1).eventtype = 'artefact_eyeblink';
S.trialdef(1).eventvalue = 'HEOG';
S.prefix = 'blink';
D_blink_epochs = spm_eeg_epochs(S); 

S = []; 
S.D = D_blink_epochs; 
S.prefix = 'ave';
D_before = spm_eeg_average(S); % average eyeblink

% spatial confounds
S = []; 
S.D = D_before; 
S.mode = 'SVD';  % Single Value Decomposition; finds the dominant spatial pattern explaining the blink topography across electrodes
S.timewin = [-50 150]; 
S.ncomp = 3; %to change with the right number of components! 1:blinks, 2:also horizontal eye movement eg; 3:also lid movement eg
D_conf = spm_eeg_spatial_confounds(S); % Kathi changed this from ...confounds_jh()

% Blink component is added to continuous dataset
S = []; 
S.D = D;  % Changed to S.D = D.ebf, from S.D = D; !!!
S.mode = 'SPMEEG';
S.conffile = D_conf;
S.prefix = 'corr1';
D_corrected = spm_eeg_spatial_confounds(S); 

% Correction: continuous EEG is projected into the subspace orthogonal to
% the blink component
S = []; 
S.D = D_corrected; 
S.mode = 'SSP'; 
S.prefix = 'corr2';
D = spm_eeg_correct_sensor_data(S); 

%display_SPM_data(D)

%% END - END EYE BLINK ON CONTINUOUS DATA


%% 3. Experimental design
% --------------------------------------------------------------------

%trial definition 
S = []; 
S.D = D; 
S.timewin = [-100 500];
S.trialdef(1).conditionlabel = 'high_intensity'; 
S.trialdef(2).conditionlabel = 'low_intensity'; 
S.trialdef(1).eventtype = 'STATUS';
S.trialdef(2).eventtype = 'STATUS';
S.trialdef(1).eventvalue = 1;
S.trialdef(2).eventvalue = 2;
S.trialdef(1).trlshift = 0; 
S.trialdef(2).trlshift = 0; 
S.reviewtrials = 0; 
S.save = 1; 
[trl, conditionlabels, S] = spm_eeg_definetrial(S); 

% new part, custom trial definition roving oddball paradigm
[new_trl, new_conditionlabels, condition_table] = classify_roving_trials(trl, conditionlabels);

% Save the condition table 
save(fullfile(D.path, 'roving_condition_table.mat'), ... 
    'condition_table');

% EPOCHING
 
%new part, epoching for roving paradigm
S = []; 
S.D = D; 
S.trl = new_trl;
S.conditionlabels = new_conditionlabels; 
S.prefix = 'roving';
D = spm_eeg_epochs(S); 


%% END EPOCHING


%% 4. Artefact detection, baseline correction
% --------------------------------------------------------------------

S = []; 
S.D = D; 
S.methods.channels = {'EEG'};
S.methods.fun = 'zscore'; 
S.methods.settings.threshold = 6; 
%S.methods.settings.excwin = 500; 
D = spm_eeg_artefact(S); 


%check the bad segments (does it make sense?)
display_SPM_data(D)


% Baseline correction
S = [];
S.D = D;
S.timewin = [-100 0];
D = spm_eeg_bc(S);

% Save the cleaned, baseline-corrected single-trial dataset
% This is the dataset from which BOTH the 20-condition and 4-condition averages will be generated.

D_single = D;


%% 5. AVERAGING
% --------------------------------------------------------------------
% Averaging is done after artifact detection and baseline correction. Trials marked as bad will be excluded from averaging.
% Output is 2 separate files, one containing all 20 conditions, the other one containing 4 conditions (SL,SH,DL,DH).

% Average 20 conditions
S = [];
S.D = D_single;
S.prefix = 'ave20';
D20 = spm_eeg_average(S);

fprintf('\n20-condition average created:\n');
fprintf('   %s\n\n', fullfile(D20.path, D20.fname));

% Create 4-condition version
% Original condition labels contain: standard_low_n2, standard_low_n3 etc.

four_conditionlabels = cell(size(new_conditionlabels));

for i = 1:length(new_conditionlabels)

    current_label = new_conditionlabels{i};

    if startsWith(current_label, 'standard_low')
        four_conditionlabels{i} = 'standard_low';

    elseif startsWith(current_label, 'standard_high')
        four_conditionlabels{i} = 'standard_high';

    elseif startsWith(current_label, 'deviant_low')
        four_conditionlabels{i} = 'deviant_low';

    elseif startsWith(current_label, 'deviant_high')
        four_conditionlabels{i} = 'deviant_high';

    else
        error('Unknown condition label: %s', current_label);
    end

end

% Create 4-condition dataset
% We now create a copy of the single-trial dataset and replace its condition assignments with the four collapsed conditions.

D4 = D_single;

% Find trials belonging to each of the 4 conditions
idx_standard_low = strcmp(four_conditionlabels, 'standard_low');
idx_standard_high = strcmp(four_conditionlabels, 'standard_high');
idx_deviant_low = strcmp(four_conditionlabels, 'deviant_low');
idx_deviant_high = strcmp(four_conditionlabels, 'deviant_high');


% Display trial counts
fprintf('\nTrial counts for 4-condition average (including bad trials):\n');

fprintf('standard_low  : %d trials\n', ...
    sum(idx_standard_low));

fprintf('standard_high : %d trials\n', ...
    sum(idx_standard_high));

fprintf('deviant_low   : %d trials\n', ...
    sum(idx_deviant_low));

fprintf('deviant_high  : %d trials\n', ...
    sum(idx_deviant_high));

% Assign the four condition labels to the individual trials!
for i = 1:D4.ntrials
    D4 = conditions(D4, i, four_conditionlabels{i});
end

% Save the modified single-trial dataset
D4.save();

% AVERAGE THE FOUR CONDITIONS
S = [];
S.D = D4;
S.prefix = 'ave4';

D4_average = spm_eeg_average(S);

% DISPLAY FINAL RESULTS

fprintf('\n');
fprintf('\n20-condition average:\n');
fprintf('   %s\n', fullfile(D20.path, D20.fname));

fprintf('\n4-condition average:\n');
fprintf('   %s\n', fullfile(D4_average.path, D4_average.fname));

fprintf('\n4 final conditions:\n');

final_conditions = condlist(D4_average);

for i = 1:length(final_conditions)

    fprintf('   %s\n', final_conditions{i});

end

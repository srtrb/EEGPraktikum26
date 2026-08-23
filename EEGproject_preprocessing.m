%% EEG PREPROCESSING PIPELINE
%{
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
%}

clear; clc; 

% Set path
project_root = ['C:\XXX'];
spm_path = ['C:\XXX'];
addpath(spm_path);
spm('defaults', 'EEG');
addpath(fullfile(spm_path, 'external', 'fieldtrip'));
cd(project_root);

%% 1. Continuous preprocessing

%convert
S = []; 
S.dataset = fullfile(project_root, '01EEG', 'raw', 'SPNCartoons_ID04.bdf');
D = spm_eeg_convert(S);

%display_SPM_data(D)

% select channels (added now for completeness)
load(fullfile(project_root, '01EEG', 'spm', 'channelselection.mat'));
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
S.freq = 0.1; 
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

load(fullfile(project_root, '01EEG', 'spm', 'avref_eog.mat'));
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

%trial definition 
S = []; 
S.D = D; 
S.timewin = [-100 500];
S.trialdef(1).conditionlabel = 'high_intensity'; 
S.trialdef(2).conditionlabel = 'low_intensity'; 
%S.trialdef(3).conditionlabel = 'block';
S.trialdef(1).eventtype = 'STATUS';
S.trialdef(2).eventtype = 'STATUS';
%S.trialdef(3).eventtype = 'STATUS';
S.trialdef(1).eventvalue = 1;
S.trialdef(2).eventvalue = 2;
%S.trialdef(3).eventvalue = 124;
S.trialdef(1).trlshift = 0; 
S.trialdef(2).trlshift = 0; 
%S.trialdef(3).trlshift = 0; 
S.reviewtrials = 0; 
S.save = 1; 
[trl, conditionlabels, S] = spm_eeg_definetrial(S); 

% new part, custom trial definition roving oddball paradigm
[new_trl, new_conditionlabels] = classify_roving_trials(trl, conditionlabels);


% EPOCHING
 
%new part, epoching for roving paradigm
S = []; 
S.D = D; 
S.trl = new_trl;
S.conditionlabels = new_conditionlabels; 
S.prefix = 'roving';
D = spm_eeg_epochs(S); 


%% END EPOCH

%% 4. Artefact detection, baseline correction, averaging

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


%average!
S = []; 
S.D = D; 
D = spm_eeg_average(S); 
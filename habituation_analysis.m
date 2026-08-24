clear; clc;

%% Load trial-level dataset
D = spm_eeg_load('barovingcorr2fMinterpolate_dfcspmeeg_SPNCartoons_ID04.mat');

nTrials = D.ntrials;
fprintf('Number of trials: %d\n\n', nTrials);

%% Conditions
% Keep original condition labels from the trial-level dataset
conditionlabels = D.conditions;

% Convert to Low / High
Condition = cell(nTrials,1);

for t = 1:nTrials
    if contains(conditionlabels{t}, 'low', 'IgnoreCase', true)
        Condition{t} = 'Low';
    elseif contains(conditionlabels{t}, 'high', 'IgnoreCase', true)
        Condition{t} = 'High';
    else
        error('Unknown condition at trial %d: %s', t, conditionlabels{t});
    end
end

Condition = categorical(Condition);

%% Repetition position
% Repetition is calculated separately within each condition.
Repeat = zeros(nTrials,1);

previous_condition = '';
counter = 0;

for t = 1:nTrials

    current_condition = Condition(t);

    if t == 1 || ~strcmp(char(current_condition), previous_condition)

        counter = 1;

    else

        counter = counter + 1;

    end

    Repeat(t) = counter;

    previous_condition = char(current_condition);

end

%% ROI definition
% Single midline ROI
ROI_channels = {'Fz','Cz','Pz'};

chan_idx = find(ismember(strtrim(D.chanlabels), ROI_channels));

if length(chan_idx) ~= 3
    error('Could not find all three ROI electrodes.');
end

fprintf('ROI: Fz + Cz + Pz\n');

for i = 1:length(chan_idx)
    fprintf('  %s -> channel %d\n', ...
        ROI_channels{i}, chan_idx(i));
end

%% Time windows
time_windows = [
    0.100 0.200
    0.200 0.350
    0.350 0.600
];

win_names = {
    'Early'
    'Attentional'
    'P300'
};

%% Main analysis

for w = 1:3

    fprintf('\n\n=============================================\n');
    fprintf('%s: %.0f-%.0f ms\n', ...
        win_names{w}, ...
        time_windows(w,1)*1000, ...
        time_windows(w,2)*1000);
    fprintf('=============================================\n');

    %% Time window indices
    erp_idx = find( ...
        D.time >= time_windows(w,1) & ...
        D.time <= time_windows(w,2));

    %% Calculate ROI amplitude for every trial
    %
    % First average over time,
    % then average Fz/Cz/Pz.
    %
    ROI_amplitude = zeros(nTrials,1);

    for t = 1:nTrials

        trial_data = D(chan_idx, erp_idx, t);

        % Mean over time
        electrode_amplitudes = mean(trial_data, 2);

        % Mean across Fz, Cz, Pz
        ROI_amplitude(t) = mean(electrode_amplitudes);

    end

    %% Create analysis table
    T = table( ...
        ROI_amplitude, ...
        Repeat, ...
        Condition, ...
        (1:nTrials)', ...
        'VariableNames', ...
        {'Amp','Repeat','Condition','Trial'});

    %% -----------------------------------------
    % LOW habituation
    % ------------------------------------------
    fprintf('\n>>> LOW habituation\n\n');

    T_low = T(T.Condition == 'Low', :);

    LME_Low = fitlme( ...
        T_low, ...
        'Amp ~ Repeat + (1|Trial)');

    disp(LME_Low.Coefficients);

    %% -----------------------------------------
    % HIGH habituation
    % ------------------------------------------
    fprintf('\n>>> HIGH habituation\n\n');

    T_high = T(T.Condition == 'High', :);

    LME_High = fitlme( ...
        T_high, ...
        'Amp ~ Repeat + (1|Trial)');

    disp(LME_High.Coefficients);

    %% -----------------------------------------
    % LOW vs HIGH comparison
    % ------------------------------------------
    fprintf('\n>>> LOW vs HIGH comparison\n\n');

    LME_Condition = fitlme( ...
        T, ...
        'Amp ~ Condition + (1|Trial)');

    disp(LME_Condition.Coefficients);

end

%% Finish
fprintf('\n\n=============================================\n');
fprintf('HABITUATION ANALYSIS COMPLETE\n');
fprintf('=============================================\n');

fprintf('Dataset: trial-level BEFORE averaging\n');
fprintf('Trials: %d\n', nTrials);
fprintf('Conditions: Low / High\n');
fprintf('ROI: Fz + Cz + Pz\n');
fprintf('Windows: 100-200, 200-350, 350-600 ms\n');
fprintf('Analysis: ROI-level only\n');
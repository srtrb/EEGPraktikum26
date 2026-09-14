%% HISTORY / TRANSITION EFFECT ANALYSIS
% SPNCartoons_ID04
%
% Main questions:
%
% 1. Does the length of the previous Low block affect
%    the ERP amplitude of the first High stimulus?
%
% 2. Does the length of the previous High block affect
%    the ERP amplitude of the first Low stimulus?
%
% Only the FIRST stimulus after an intensity transition is analysed.
%
% ROI: Fz + FCz + Cz + C3 + C4
% ERP window: 0-150 ms

clearvars;
clc;

%% DATASET

D = spm_eeg_load( ...
    'barovingcorr2fMinterpolate_dfcspmeeg_SPNCartoons_ID04.mat');

nTrials = D.ntrials;

fprintf('Number of trials: %d\n\n', nTrials);

%% ELECTRODES

roi_labels = {'Fz','FCz','Cz','C3','C4'};

roi_channels = zeros(1, numel(roi_labels));

for i = 1:numel(roi_labels)

    roi_channels(i) = D.indchannel(roi_labels{i});

end

fprintf('ROI: Fz + FCz + Cz + C3 + C4\n');

for i = 1:numel(roi_labels)

    fprintf('  %s -> channel %d\n', ...
        roi_labels{i}, roi_channels(i));

end

%% CONDITIONS -> LOW / HIGH

conditions = D.conditions;

Condition = strings(nTrials,1);

for t = 1:nTrials

    cond = lower(string(conditions{t}));

    if contains(cond,'low')

        Condition(t) = "Low";

    elseif contains(cond,'high')

        Condition(t) = "High";

    else

        Condition(t) = missing;

    end

end

%% IDENTIFY TRANSITIONS AND PREVIOUS BLOCK LENGTH

PrevCondition = strings(nTrials,1);

PreviousCount = nan(nTrials,1);

IsTransition = false(nTrials,1);

current_condition = "";
current_run_length = 0;

for t = 1:nTrials

    if ismissing(Condition(t))

        current_condition = "";
        current_run_length = 0;

        continue;

    end

    % First valid trial
    if current_condition == ""

        current_condition = Condition(t);
        current_run_length = 1;

    else

        % Same intensity as previous trial
        if Condition(t) == current_condition

            current_run_length = current_run_length + 1;

        % Intensity changed
        else

            % Current trial is the first trial of the new block
            PrevCondition(t) = current_condition;

            PreviousCount(t) = current_run_length;

            IsTransition(t) = true;

            % Start counting the new block
            current_condition = Condition(t);
            current_run_length = 1;

        end

    end

end

%% TRANSITION SUMMARY

fprintf('\n');
fprintf('=============================================\n');
fprintf('TRANSITION SUMMARY\n');
fprintf('=============================================\n');

idx_LowHigh = IsTransition & ...
              PrevCondition == "Low" & ...
              Condition == "High";

idx_HighLow = IsTransition & ...
              PrevCondition == "High" & ...
              Condition == "Low";

fprintf('Low -> High transitions: %d\n', ...
    sum(idx_LowHigh));

fprintf('High -> Low transitions: %d\n', ...
    sum(idx_HighLow));

%% ERP DATA

time = D.time;

time_ms = time * 1000;

% Select ROI channels
X = D(roi_channels,:,:);

% Average across ROI channels
X = squeeze(mean(X,1));

% X = time x trials

fprintf('\nERP time range: %.0f to %.0f ms\n', ...
    time_ms(1), time_ms(end));

%% ERP WINDOW

idx_window = time_ms >= 0 & time_ms <= 150;

fprintf('ERP analysis window: 0 to 150 ms\n');

%% =========================================================
% LOW -> HIGH
% ==========================================================

fprintf('\n');
fprintf('=============================================\n');
fprintf('LOW -> HIGH\n');
fprintf('=============================================\n');

% First High trial after a Low block
idx = idx_LowHigh;

Count_LH = PreviousCount(idx);

% ERP amplitude of each individual transition trial
ERP_LH = zeros(sum(idx),1);

trial_indices_LH = find(idx);

for i = 1:length(trial_indices_LH)

    trial = trial_indices_LH(i);

    ERP_LH(i) = mean(X(idx_window,trial));

end

% Remove groups with fewer than 3 trials
unique_counts_LH = unique(Count_LH);

valid_counts_LH = [];

for c = 1:length(unique_counts_LH)

    n = sum(Count_LH == unique_counts_LH(c));

    if n >= 3
        valid_counts_LH(end+1) = unique_counts_LH(c);
    end

end

keep_LH = ismember(Count_LH, valid_counts_LH);

Count_LH = Count_LH(keep_LH);
ERP_LH = ERP_LH(keep_LH);

%% TABLE: LOW -> HIGH

T_LH = table( ...
    Count_LH, ...
    ERP_LH, ...
    'VariableNames', ...
    {'PreviousCount','ERP_Amplitude'});

disp(T_LH);

%% REGRESSION: LOW -> HIGH

mdl_LH = fitlm(Count_LH, ERP_LH);

fprintf('\nRegression: Low -> High\n');

disp(mdl_LH);

fprintf('\nRegression coefficients:\n');

disp(mdl_LH.Coefficients);

%% =========================================================
% HIGH -> LOW
% ==========================================================

fprintf('\n');
fprintf('=============================================\n');
fprintf('HIGH -> LOW\n');
fprintf('=============================================\n');

% First Low trial after a High block
idx = idx_HighLow;

Count_HL = PreviousCount(idx);

% ERP amplitude of each individual transition trial
ERP_HL = zeros(sum(idx),1);

trial_indices_HL = find(idx);

for i = 1:length(trial_indices_HL)

    trial = trial_indices_HL(i);

    ERP_HL(i) = mean(X(idx_window,trial));

end

% Remove groups with fewer than 3 trials
unique_counts_HL = unique(Count_HL);

valid_counts_HL = [];

for c = 1:length(unique_counts_HL)

    n = sum(Count_HL == unique_counts_HL(c));

    if n >= 3
        valid_counts_HL(end+1) = unique_counts_HL(c);
    end

end

keep_HL = ismember(Count_HL, valid_counts_HL);

Count_HL = Count_HL(keep_HL);
ERP_HL = ERP_HL(keep_HL);

%% TABLE: HIGH -> LOW

T_HL = table( ...
    Count_HL, ...
    ERP_HL, ...
    'VariableNames', ...
    {'PreviousCount','ERP_Amplitude'});

disp(T_HL);

%% REGRESSION: HIGH -> LOW

mdl_HL = fitlm(Count_HL, ERP_HL);

fprintf('\nRegression: High -> Low\n');

disp(mdl_HL);

fprintf('\nRegression coefficients:\n');

disp(mdl_HL.Coefficients);

%% =========================================================
% SUMMARY
% ==========================================================

fprintf('\n');
fprintf('=============================================\n');
fprintf('SUMMARY\n');
fprintf('=============================================\n');

Summary = table( ...
    ["Low -> High"; "High -> Low"], ...
    [mdl_LH.Coefficients.Estimate(2); ...
     mdl_HL.Coefficients.Estimate(2)], ...
    [mdl_LH.Coefficients.SE(2); ...
     mdl_HL.Coefficients.SE(2)], ...
    [mdl_LH.Coefficients.pValue(2); ...
     mdl_HL.Coefficients.pValue(2)], ...
    'VariableNames', ...
    {'Transition','Beta','SE','pValue'});

disp(Summary);

fprintf('\n');
fprintf('=============================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('=============================================\n');
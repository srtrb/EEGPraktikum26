%% HABITUATION ANALYSIS
% Trial-level, before averaging
% ROI: Fz + Cz + Pz

clearvars -except D
clc

%% ---------------------------------------------------------
% DATASET
%% ---------------------------------------------------------

D = spm_eeg_load('barovingcorr2fMinterpolate_dfcspmeeg_SPNCartoons_ID04.mat');

nTrials = D.ntrials;

fprintf('Number of trials: %d\n\n', nTrials);

%% ---------------------------------------------------------
% ELECTRODES
%% ---------------------------------------------------------

roi_labels = {'Fz','Cz','Pz'};
roi_channels = zeros(1,3);

for i = 1:3
    roi_channels(i) = D.indchannel(roi_labels{i});
end

fprintf('ROI: Fz + Cz + Pz\n');
for i = 1:3
    fprintf('  %s -> channel %d\n', ...
        roi_labels{i}, roi_channels(i));
end

%% ---------------------------------------------------------
% CONDITIONS -> LOW / HIGH
%% ---------------------------------------------------------

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

%% ---------------------------------------------------------
% REPETITION POSITION
% Number of consecutive previous trials with same condition
%% ---------------------------------------------------------

Repeat = nan(nTrials,1);

current_repeat = 0;
previous_condition = "";

for t = 1:nTrials

    if ismissing(Condition(t))
        current_repeat = 0;
        previous_condition = "";
        continue
    end

    if Condition(t) == previous_condition
        current_repeat = current_repeat + 1;
    else
        current_repeat = 1;
    end

    Repeat(t) = current_repeat;

    previous_condition = Condition(t);
end

%% ---------------------------------------------------------
% TIME WINDOWS
%% ---------------------------------------------------------

windows = {
    'Early',       [100 200];
    'Attentional', [200 350];
    'P300',        [350 600]
};

%% ---------------------------------------------------------
% ANALYSIS
%% ---------------------------------------------------------

for w = 1:size(windows,1)

    win_name = windows{w,1};
    timewin  = windows{w,2};

    fprintf('\n');
    fprintf('=============================================\n');
    fprintf('%s: %d-%d ms\n', ...
        win_name, timewin(1), timewin(2));
    fprintf('=============================================\n');

    % Find samples using D.time
    time_idx = find(D.time >= timewin(1)/1000 & ...
                    D.time <= timewin(2)/1000);

    if isempty(time_idx)
        fprintf('No samples found for this window.\n');
        continue
    end

    %% -----------------------------------------------------
    % Extract amplitude
    %% -----------------------------------------------------

    Y = nan(nTrials,1);

    for t = 1:nTrials

        if ismissing(Condition(t))
            continue
        end

        % ROI x time
        x = D(roi_channels,time_idx,t);

        % Mean across electrodes and time
        Y(t) = mean(x(:),'omitnan');
    end

    %% -----------------------------------------------------
    % LOW
    %% -----------------------------------------------------

    idx = Condition == "Low" & ~isnan(Y) & ~isnan(Repeat);

    T_low = table( ...
        Y(idx), ...
        Repeat(idx), ...
        'VariableNames', {'Amplitude','Repeat'});

    fprintf('\n>>> LOW habituation\n\n');

    if numel(unique(T_low.Repeat)) < 2

        fprintf('Not enough repetition variation.\n');

    else

        lme_low = fitlme( ...
            T_low, ...
            'Amplitude ~ Repeat');

        disp(lme_low.Coefficients);
    end

    %% -----------------------------------------------------
    % HIGH
    %% -----------------------------------------------------

    idx = Condition == "High" & ~isnan(Y) & ~isnan(Repeat);

    T_high = table( ...
        Y(idx), ...
        Repeat(idx), ...
        'VariableNames', {'Amplitude','Repeat'});

    fprintf('\n>>> HIGH habituation\n\n');

    if numel(unique(T_high.Repeat)) < 2

        fprintf('Not enough repetition variation.\n');

    else

        lme_high = fitlme( ...
            T_high, ...
            'Amplitude ~ Repeat');

        disp(lme_high.Coefficients);
    end

    %% -----------------------------------------------------
    % LOW vs HIGH
    %% -----------------------------------------------------

    idx = ~isnan(Y) & ~ismissing(Condition);

    T = table( ...
        Y(idx), ...
        categorical(Condition(idx)), ...
        'VariableNames', {'Amplitude','Condition'});

    fprintf('\n>>> LOW vs HIGH comparison\n\n');

    lme_condition = fitlme( ...
        T, ...
        'Amplitude ~ Condition');

    disp(lme_condition.Coefficients);

end

fprintf('\n');
fprintf('=============================================\n');
fprintf('HABITUATION ANALYSIS COMPLETE\n');
fprintf('=============================================\n');

fprintf('Dataset: trial-level BEFORE averaging\n');
fprintf('Trials: %d\n', nTrials);
fprintf('Conditions: Low / High\n');
fprintf('ROI: Fz + Cz + Pz\n');
fprintf('Time windows: 100-200, 200-350, 350-600 ms\n');
fprintf('Analysis: ALL repetitions\n');
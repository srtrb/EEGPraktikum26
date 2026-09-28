%% ADAPTATION ANALYSIS
% ID04
%
% ROI: C2, C4, CP2, CP4
% Time window: 0-200 ms
%
% Low -> High:
% ERP amplitude of first HIGH stimulus
% as a function of preceding LOW rove length
%
% High -> Low:
% ERP amplitude of first LOW stimulus
% as a function of preceding HIGH rove length


clear;
clc;


%% ------------------------------------------------------------------------
% 1. LOAD DATA
% -------------------------------------------------------------------------

data_path = 'C:\Users\MCNB01\Documents\MATLAB\data\EEGPraktikum26\Group4';

eeg_file = fullfile(data_path, ...
    'barovingcorr2fMinterpolate_dfcspmeeg_SPNCartoons_ID04.mat');

trialdef_file = fullfile(data_path, 'trialdef.mat');

D = spm_eeg_load(eeg_file);

load(trialdef_file, 'conditionlabels');


%% ------------------------------------------------------------------------
% 2. RECONSTRUCT ROVES
% -------------------------------------------------------------------------

% Make condition labels a column vector
conditionlabels = conditionlabels(:);

% Low = 1
% High = 2

values = zeros(numel(conditionlabels),1);

values(strcmp(conditionlabels,'low_intensity')) = 1;
values(strcmp(conditionlabels,'high_intensity')) = 2;

if any(values == 0)
    error('Unexpected condition label found.');
end

rove_start = [ ...
    1
    find(diff(values) ~= 0) + 1
];

rove_start = rove_start(:);

rove_end = [ ...
    rove_start(2:end) - 1
    numel(values)
];

rove_end = rove_end(:);

n_roves = numel(rove_start);


%% ------------------------------------------------------------------------
% 3. RECONSTRUCT RETAINED TRIALS
% -------------------------------------------------------------------------

keep_idx    = [];
condition   = {};
intensity   = {};
rove_length = [];
transition  = {};

for r = 1:n_roves

    start_idx = rove_start(r);
    end_idx   = rove_end(r);

    current_rove_length = end_idx - start_idx + 1;

    current_label = conditionlabels{start_idx};

    if strcmp(current_label,'low_intensity')

        current_intensity = 'low';

    elseif strcmp(current_label,'high_intensity')

        current_intensity = 'high';

    else

        error('Unexpected condition label.');

    end


    % ---------------------------------------------------------------------
    % STANDARD
    % Last trial of current rove
    % ---------------------------------------------------------------------

    if current_rove_length >= 2

        keep_idx(end+1,1) = end_idx;

        condition{end+1,1} = 'standard';
        intensity{end+1,1} = current_intensity;
        rove_length(end+1,1) = current_rove_length;
        transition{end+1,1} = 'standard';

    end


    % ---------------------------------------------------------------------
    % DEVIANT
    % First trial of next rove
    % ---------------------------------------------------------------------

    if r < n_roves && current_rove_length >= 2

        deviant_idx = rove_start(r+1);

        deviant_label = conditionlabels{deviant_idx};

        if strcmp(deviant_label,'low_intensity')

            deviant_intensity = 'low';
            deviant_transition = 'High_to_Low';

        elseif strcmp(deviant_label,'high_intensity')

            deviant_intensity = 'high';
            deviant_transition = 'Low_to_High';

        else

            error('Unexpected deviant condition label.');

        end


        keep_idx(end+1,1) = deviant_idx;

        condition{end+1,1} = 'deviant';
        intensity{end+1,1} = deviant_intensity;

        % Deviant is associated with the preceding rove length
        rove_length(end+1,1) = current_rove_length;

        transition{end+1,1} = deviant_transition;

    end

end


%% ------------------------------------------------------------------------
% 4. SORT BY ORIGINAL TRIAL ORDER
% -------------------------------------------------------------------------

[keep_idx, order] = sort(keep_idx);

condition   = condition(order);
intensity   = intensity(order);
rove_length = rove_length(order);
transition  = transition(order);


%% ------------------------------------------------------------------------
% 5. CHECK RECONSTRUCTION
% -------------------------------------------------------------------------

if numel(keep_idx) ~= D.ntrials

    error(['Reconstructed %d trials, but EEG contains %d trials.'], ...
        numel(keep_idx), D.ntrials);

end


fprintf('Retained trials: %d\n', numel(keep_idx));


%% ------------------------------------------------------------------------
% 6. def ROI
% -------------------------------------------------------------------------

roi_labels = {'C2', 'C4', 'CP2', 'CP4'};

roi_channels = zeros(1,numel(roi_labels));

for i = 1:numel(roi_labels)

    roi_channels(i) = indchannel(D,roi_labels{i});

end


%% ------------------------------------------------------------------------
% 7. TIME WINDOW: 0-200 ms
% -------------------------------------------------------------------------

time_ms = D.time(:) * 1000;

sample_idx = find( ...
    time_ms >= 0 & ...
    time_ms <= 150);

fprintf('Time window: %.0f-%.0f ms\n', ...
    time_ms(sample_idx(1)), ...
    time_ms(sample_idx(end)));

fprintf('Number of samples: %d\n', ...
    numel(sample_idx));


%% ------------------------------------------------------------------------
% 8. EXTRACT ERP AMPLITUDE
% -------------------------------------------------------------------------

% Channels × time × trials

Y = double(D( ...
    roi_channels, ...
    sample_idx, ...
    1:D.ntrials));


% Average across ROI channels

Y = squeeze(mean(Y,1));


% Average across 0-200 ms

amplitude = squeeze(mean(Y,1));

amplitude = amplitude(:);


%% ------------------------------------------------------------------------
% 9. SELECT DEVIANT TRIALS
% -------------------------------------------------------------------------

is_deviant = strcmp(condition,'deviant');

amplitude = amplitude(is_deviant);
rove_length = rove_length(is_deviant);
transition = transition(is_deviant);


%% ------------------------------------------------------------------------
% 10. LOW -> HIGH
% -------------------------------------------------------------------------

idx_LH = strcmp(transition,'Low_to_High');

X_LH = rove_length(idx_LH);
Y_LH = amplitude(idx_LH);

fprintf('Low -> High: %d trials\n',numel(Y_LH));


%% ------------------------------------------------------------------------
% 11. HIGH -> LOW
% -------------------------------------------------------------------------

idx_HL = strcmp(transition,'High_to_Low');

X_HL = rove_length(idx_HL);
Y_HL = amplitude(idx_HL);

fprintf('High -> Low: %d trials\n',numel(Y_HL));


%% ------------------------------------------------------------------------
% 12. LINEAR REGRESSION
% -------------------------------------------------------------------------

mdl_LH = fitlm(X_LH,Y_LH);

mdl_HL = fitlm(X_HL,Y_HL);


%% ------------------------------------------------------------------------
% 13. EXTRACT RESULTS
% -------------------------------------------------------------------------

beta_LH = mdl_LH.Coefficients.Estimate(2);
SE_LH   = mdl_LH.Coefficients.SE(2);
t_LH    = mdl_LH.Coefficients.tStat(2);
p_LH    = mdl_LH.Coefficients.pValue(2);
R2_LH   = mdl_LH.Rsquared.Ordinary;


beta_HL = mdl_HL.Coefficients.Estimate(2);
SE_HL   = mdl_HL.Coefficients.SE(2);
t_HL    = mdl_HL.Coefficients.tStat(2);
p_HL    = mdl_HL.Coefficients.pValue(2);
R2_HL   = mdl_HL.Rsquared.Ordinary;


%% ------------------------------------------------------------------------
% 14. RESULTS
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf('HABITUATION REGRESSION RESULTS\n');
fprintf('============================================================\n');

fprintf('\nLOW -> HIGH\n');
fprintf('N    = %d\n',numel(Y_LH));
fprintf('Beta = %.6f\n',beta_LH);
fprintf('SE   = %.6f\n',SE_LH);
fprintf('t    = %.6f\n',t_LH);
fprintf('p    = %.6g\n',p_LH);
fprintf('R^2  = %.6f\n',R2_LH);

fprintf('\nHIGH -> LOW\n');
fprintf('N    = %d\n',numel(Y_HL));
fprintf('Beta = %.6f\n',beta_HL);
fprintf('SE   = %.6f\n',SE_HL);
fprintf('t    = %.6f\n',t_HL);
fprintf('p    = %.6g\n',p_HL);
fprintf('R^2  = %.6f\n',R2_HL);

fprintf('\n============================================================\n');
fprintf('ROI: C2, C4, CP2, CP4\n');
fprintf('Time window: 0-200 ms\n');
fprintf('============================================================\n');
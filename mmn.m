%% MMN ANALYSIS
% SPNCartoons_ID04
%
% Dataset: 4-condition averaged ERP
% ROI: C2 + CP2 + CP4 + FC6
% MMN: Deviant - Standard
% Plot: 0-200 ms
% Analysis window: 50-200 ms

clearvars;
clc;

%% ------------------------------------------------------------------------
% DATASET
% -------------------------------------------------------------------------

D = spm_eeg_load( ...
    'ave4barovingcorr2fMinterpolate_dfcspmeeg_SPNCartoons_ID04.mat');


%% ------------------------------------------------------------------------
% ROI
% -------------------------------------------------------------------------
roi_labels = {'C2','CP2','CP4','FC6'};
% roi_labels = {'Fz','FCz','Cz'};

roi_channels = zeros(1, numel(roi_labels));

for i = 1:numel(roi_labels)

    roi_channels(i) = D.indchannel(roi_labels{i});

end

fprintf('\nROI channels:\n');
disp(roi_labels);

fprintf('Channel indices:\n');
disp(roi_channels);


%% ------------------------------------------------------------------------
% CONDITIONS
% -------------------------------------------------------------------------

conditions = string(D.conditions);
conditions = conditions(:);

idx_standard_low  = find(conditions == "standard_low");
idx_standard_high = find(conditions == "standard_high");
idx_deviant_low   = find(conditions == "deviant_low");
idx_deviant_high  = find(conditions == "deviant_high");

fprintf('\nCondition indices:\n');

fprintf('standard_low  : %d\n', idx_standard_low);
fprintf('standard_high : %d\n', idx_standard_high);
fprintf('deviant_low   : %d\n', idx_deviant_low);
fprintf('deviant_high  : %d\n', idx_deviant_high);



%% ------------------------------------------------------------------------
% TIME
% -------------------------------------------------------------------------

time = D.time;
time_ms = time * 1000;

fprintf('\nTime range:\n');
fprintf('%8.2f to %8.2f ms\n', ...
    min(time_ms), max(time_ms));


%% ------------------------------------------------------------------------
% DATA
% -------------------------------------------------------------------------

% Extract ROI data

X = D(roi_channels,:,:);

fprintf('\nRaw extracted data dimensions:\n');
disp(size(X));


%% ------------------------------------------------------------------------
% AVERAGE ACROSS ROI CHANNELS
% -------------------------------------------------------------------------

% Average across C2, CP2, CP4 and FC6

X_roi = mean(X,1);

fprintf('After ROI averaging:\n');
disp(size(X_roi));


%% ------------------------------------------------------------------------
% EXTRACT ERP FOR EACH CONDITION
% -------------------------------------------------------------------------

ERP_StandardLow = squeeze( ...
    X_roi(:,:,idx_standard_low));

ERP_StandardHigh = squeeze( ...
    X_roi(:,:,idx_standard_high));

ERP_DeviantLow = squeeze( ...
    X_roi(:,:,idx_deviant_low));

ERP_DeviantHigh = squeeze( ...
    X_roi(:,:,idx_deviant_high));


% Force column vectors

ERP_StandardLow  = ERP_StandardLow(:);
ERP_StandardHigh = ERP_StandardHigh(:);
ERP_DeviantLow   = ERP_DeviantLow(:);
ERP_DeviantHigh  = ERP_DeviantHigh(:);


%% ------------------------------------------------------------------------
% CHECK DATA LENGTHS
% -------------------------------------------------------------------------

fprintf('\nERP vector lengths:\n');

fprintf('Standard Low  : %d\n', numel(ERP_StandardLow));
fprintf('Standard High : %d\n', numel(ERP_StandardHigh));
fprintf('Deviant Low   : %d\n', numel(ERP_DeviantLow));
fprintf('Deviant High  : %d\n', numel(ERP_DeviantHigh));
fprintf('Time points   : %d\n', numel(time_ms));


%% ------------------------------------------------------------------------
% CHECK FOR NaN / INF
% -------------------------------------------------------------------------

fprintf('\nNaN / Inf check:\n');

fprintf('Standard Low  : %d NaN, %d Inf\n', ...
    sum(isnan(ERP_StandardLow)), ...
    sum(isinf(ERP_StandardLow)));

fprintf('Standard High : %d NaN, %d Inf\n', ...
    sum(isnan(ERP_StandardHigh)), ...
    sum(isinf(ERP_StandardHigh)));

fprintf('Deviant Low   : %d NaN, %d Inf\n', ...
    sum(isnan(ERP_DeviantLow)), ...
    sum(isinf(ERP_DeviantLow)));

fprintf('Deviant High  : %d NaN, %d Inf\n', ...
    sum(isnan(ERP_DeviantHigh)), ...
    sum(isinf(ERP_DeviantHigh)));


%% ------------------------------------------------------------------------
% OVERALL STANDARD / DEVIANT
% -------------------------------------------------------------------------

% Average across Low and High intensity

ERP_StandardOverall = ...
    (ERP_StandardLow + ERP_StandardHigh) / 2;

ERP_DeviantOverall = ...
    (ERP_DeviantLow + ERP_DeviantHigh) / 2;


%% ------------------------------------------------------------------------
% MMN
% -------------------------------------------------------------------------

% Overall MMN = Deviant - Standard

MMN_Overall = ...
    ERP_DeviantOverall - ERP_StandardOverall;


%% ------------------------------------------------------------------------
% CHECK NUMERICAL RANGE
% -------------------------------------------------------------------------

fprintf('\nERP amplitude range:\n');

fprintf('Standard: %8.4f to %8.4f uV\n', ...
    min(ERP_StandardOverall), ...
    max(ERP_StandardOverall));

fprintf('Deviant : %8.4f to %8.4f uV\n', ...
    min(ERP_DeviantOverall), ...
    max(ERP_DeviantOverall));

fprintf('MMN     : %8.4f to %8.4f uV\n', ...
    min(MMN_Overall), ...
    max(MMN_Overall));


%% ------------------------------------------------------------------------
% PLOT: STANDARD / DEVIANT / MMN
% -------------------------------------------------------------------------

figure('Color','w');

% Standard
plot(time_ms, ERP_StandardOverall, ...
    'Color', [0 0.4470 0.7410], ...
    'LineWidth', 2);

hold on;

% Deviant
plot(time_ms, ERP_DeviantOverall, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 2);

% MMN
plot(time_ms, MMN_Overall, ...
    'Color', [0.4940 0.1840 0.5560], ...
    'LineWidth', 2);

% Zero line
yline(0,'k--');

grid on;
box on;

% Display only 0-200 ms
xlim([0 200]);

% EEG convention
set(gca,'YDir','reverse');

xlabel('Time (ms)');
ylabel('Amplitude (\muV)');

title('Mismatch Negativity');

legend({ ...
    'Standard', ...
    'Deviant', ...
    'MMN'}, ...
    'Location','best');


%% ------------------------------------------------------------------------
% 50-200 ms ANALYSIS
% -------------------------------------------------------------------------

idx_50_200 = ...
    time_ms >= 50 & time_ms <= 200;

standard_50_200 = ...
    mean(ERP_StandardOverall(idx_50_200));

deviant_50_200 = ...
    mean(ERP_DeviantOverall(idx_50_200));

mmn_50_200 = ...
    mean(MMN_Overall(idx_50_200));


%% ------------------------------------------------------------------------
% RESULTS
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('===============================================\n');
fprintf('OVERALL MMN: 50-200 ms\n');
fprintf('===============================================\n');

fprintf('\nROI:\n');
fprintf('C2 + CP2 + CP4 + FC6\n');

fprintf('\nTime window:\n');
fprintf('50-200 ms\n');

fprintf('\nStandard : %8.4f uV\n', ...
    standard_50_200);

fprintf('Deviant  : %8.4f uV\n', ...
    deviant_50_200);

fprintf('MMN      : %8.4f uV\n', ...
    mmn_50_200);

fprintf('\n===============================================\n');
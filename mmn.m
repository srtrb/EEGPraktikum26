%% MMN ANALYSIS
% SPNCartoons_ID04
% Analysis window: 50-200 ms

clearvars;
clc;

%% DATASET

D = spm_eeg_load( ...
    'barovingcorr2fMinterpolate_dfcspmeeg_SPNCartoons_ID04.mat');

%% ROI

% roi_labels = {'Fz,'FCz','Cz','C3','C4'};
roi_labels = {'C2','CP2','CP4','FC6'};

roi_channels = zeros(1, numel(roi_labels));

for i = 1:numel(roi_labels)
    roi_channels(i) = D.indchannel(roi_labels{i});
end

%% CONDITIONS

conditions = string(D.conditions);
conditions = conditions(:);

isStandardLow  = conditions == "standard_low";
isStandardHigh = conditions == "standard_high";

isDeviantLow   = conditions == "deviant_low";
isDeviantHigh  = conditions == "deviant_high";

%% CHECK CONDITION COUNTS

fprintf('\nCondition counts:\n');

fprintf('Standard Low  : %d\n', sum(isStandardLow));
fprintf('Standard High : %d\n', sum(isStandardHigh));
fprintf('Deviant Low   : %d\n', sum(isDeviantLow));
fprintf('Deviant High  : %d\n', sum(isDeviantHigh));

%% DATA

time = D.time;
time_ms = time * 1000;

X = D(roi_channels,:,:);

% Average across ROI channels
X_roi = squeeze(mean(X,1));

%% ERP BY INTENSITY

ERP_StandardLow = ...
    mean(X_roi(:,isStandardLow),2);

ERP_StandardHigh = ...
    mean(X_roi(:,isStandardHigh),2);

ERP_DeviantLow = ...
    mean(X_roi(:,isDeviantLow),2);

ERP_DeviantHigh = ...
    mean(X_roi(:,isDeviantHigh),2);

%% STANDARD AND DEVIANT ERP

% Equal weighting of Low and High intensity

ERP_StandardOverall = ...
    (ERP_StandardLow + ERP_StandardHigh) / 2;

ERP_DeviantOverall = ...
    (ERP_DeviantLow + ERP_DeviantHigh) / 2;

%% OVERALL MMN

MMN_Overall = ...
    ERP_DeviantOverall - ERP_StandardOverall;

%% PLOT: STANDARD, DEVIANT AND MMN

figure('Color','w');

plot(time_ms, ERP_StandardOverall, ...
    'k', 'LineWidth',2);

hold on;

plot(time_ms, ERP_DeviantOverall, ...
    'r', 'LineWidth',2);

plot(time_ms, MMN_Overall, ...
    'b', 'LineWidth',2);

yline(0, 'k--');

hold off;

grid on;
box on;

xlim([50 200]);

set(gca,'YDir','reverse');

xlabel('Time (ms)');
ylabel('Amplitude (\muV)');

title('Mismatch Negativity');

legend({ ...
    'Standard', ...
    'Deviant', ...
    'MMN'}, ...
    'Location','best');

%% 50 to 200 ms ANALYSIS

idx_50_200 = time_ms >= 50 & time_ms <= 200;

%% MEAN ERP AMPLITUDES

standard_50_200 = ...
    mean(ERP_StandardOverall(idx_50_200));

deviant_50_200 = ...
    mean(ERP_DeviantOverall(idx_50_200));

%% OVERALL MMN

mmn_50_200 = ...
    mean(MMN_Overall(idx_50_200));

%% RESULTS

fprintf('\n');
fprintf('===============================================\n');
fprintf('OVERALL MMN: 50-200 ms\n');
fprintf('===============================================\n');

fprintf('\nStandard : %8.4f uV\n', ...
    standard_50_200);

fprintf('Deviant  : %8.4f uV\n', ...
    deviant_50_200);

fprintf('MMN      : %8.4f uV\n', ...
    mmn_50_200);

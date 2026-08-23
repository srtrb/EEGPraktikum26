function [new_trl, new_conditionlabels] = classify_roving_trials(trl, conditionlabels)
% CLASSIFY_ROVING_TRIALS
%
% Classifies trials from a roving oddball paradigm into:
%   - standard_low
%   - standard_high
%   - deviant_low
%   - deviant_high
%
% Only the LAST standard before each intensity change is kept.
%
% INPUTS
%   trl               Trial matrix returned by spm_eeg_definetrial
%   conditionlabels   Cell array returned by spm_eeg_definetrial
%
% OUTPUTS
%   new_trl
%   new_conditionlabels

%% Convert condition labels to numeric values
values = zeros(length(conditionlabels),1);

values(strcmp(conditionlabels,'low_intensity'))  = 1;
values(strcmp(conditionlabels,'high_intensity')) = 2;

%% Find deviants and last standards
n = length(values);

isDeviant = false(n,1);
isLastStandard = false(n,1);

% Deviant = first trial after a change
isDeviant(2:end) = diff(values) ~= 0;

% Last standard = trial immediately before a change
isLastStandard(1:end-1) = diff(values) ~= 0;

%% Keep only deviants and last standards
keep = isDeviant | isLastStandard;

new_trl = trl(keep,:);

%% Create new condition labels
new_conditionlabels = cell(sum(keep),1);

k = 1;

for i = 1:n

    if isDeviant(i)

        if values(i)==1
            new_conditionlabels{k} = 'deviant_low';
        else
            new_conditionlabels{k} = 'deviant_high';
        end
        k = k + 1;

    elseif isLastStandard(i)

        if values(i)==1
            new_conditionlabels{k} = 'standard_low';
        else
            new_conditionlabels{k} = 'standard_high';
        end
        k = k + 1;

    end

end

end
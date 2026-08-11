%% WP4 PUPILLOMETRY PRELIMINARY ANALYSIS
% Standalone script: it does not require another script or workspace.
%
% Method follows the WP3 pupillometry logic:
%   1) invalid left/right pupil samples are removed;
%   2) the available eyes are averaged at each timestamp;
%   3) eye samples are matched many-to-one to behavioural trials by Trial ID;
%   4) time is expressed relative to tactile onset (Time_vibration);
%   5) pupil diameter is baseline-corrected using -0.4 to 0 s;
%   6) mean pupil change is calculated from 0 to 0.5 s after tactile onset;
%   7) only central trials (Monster == 0) are analysed.
%
% The behavioural CSV is used only as a trial/event log. No behavioural
% outcome (e.g., RT or performance) is analysed in this script.
%
% WP4 adaptation:
%   Baseline and Post-training are compared within participant using exact
%   sign-flip tests.
%   Five visuo-tactile distances are analysed, together with the tactile-only
%   control (distance 0). The visual-only control (distance -1) is not part
%   of the tactile-onset analysis because it has no tactile event.

clear; close all; clc;
dbstop if error
rng(2026, 'twister');

%% ========================= USER SETTINGS ==============================

% If eye and behavioural files are in the same folder, keep both paths equal.

% If the two file types are stored in separate folders with identical names,
% point each path to the corresponding folder. The script will also accept
% behavioural files named 02AE_pre.csv in a separate behavioural folder.
eyeDataFolder       = '';
behaviourDataFolder = '';
outputFolder        = fullfile(eyeDataFolder, 'WP4_Pupillometry_Outputs');

% Participant 01AM is not listed because this participant was excluded by
% the previously completed baseline PPS quality-control procedure.
participantIDs = {''};
phaseNames     = {'Baseline','Post-training'};
phaseTokens    = {'pre','post'};

% WP3 preprocessing parameters
baselineWindow = [-0.4 0];
responseWindow = [0 0.5];
minimumValidFraction = 0.50;
minimumValidTrialsPerCondition = 8; % at least 50% of the 16 planned trials

ppsDistances = [0.25 0.75 1.25 1.75 2.25];

% Tactile-aligned conditions. Distance 0 is tactile-only.
analysisDistances = [0 0.25 0.75 1.25 1.75 2.25];
conditionLabels = {'Tactile only','0.25 m','0.75 m','1.25 m','1.75 m','2.25 m'};
mainDistanceIndices = 2:6;

% Common time axis is used only for plots, not for the scalar analysis.
plotSamplingRate = 90;
commonTime = baselineWindow(1):(1/plotSamplingRate):responseWindow(2);

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ====================== INCLUDED PARTICIPANTS ========================

% The participant list already reflects the exclusions determined by the
% separate PPS behavioural quality-control analysis. This pupil script does
% not repeat or perform any behavioural analysis.
includedParticipantIDs = participantIDs;
nParticipants = numel(includedParticipantIDs);

fprintf('\nParticipants entering the pupil analysis: %d\n', nParticipants);
disp(includedParticipantIDs(:));

if nParticipants < 2
    error('Fewer than two participants remain after PPS baseline QC.');
end

%% ================== TRIAL-LEVEL PUPIL PROCESSING =====================

nPhases = numel(phaseNames);
nConditions = numel(analysisDistances);
nTime = numel(commonTime);

% Participant-level time courses: participant x phase x condition x time.
participantTimeCourse = nan(nParticipants, nPhases, nConditions, nTime);

trialRows = cell(0, 12);

for s = 1:nParticipants
    participant = includedParticipantIDs{s};

    for p = 1:nPhases
        stem = [participant '_' phaseTokens{p}];
        eyeFile = local_resolveFile(eyeDataFolder, stem, 'eye');
        behaviourFile = local_resolveFile(behaviourDataFolder, stem, 'behaviour');

        E = local_readAndCleanNames(eyeFile);
        B = local_readAndCleanNames(behaviourFile);

        eyeTrial = local_getNumericColumn(E, 'Trial');
        eyeTimestamp = local_getNumericColumn(E, 'Timestamp');
        pupilLeft = local_getNumericColumn(E, 'Pupil_diameter_L');
        pupilRight = local_getNumericColumn(E, 'Pupil_diameter_R');

        pupilLeft(~isfinite(pupilLeft) | pupilLeft <= 0) = NaN;
        pupilRight(~isfinite(pupilRight) | pupilRight <= 0) = NaN;
        pupilAverage = local_rowMeanTwoEyes(pupilLeft, pupilRight);

        behaviouralTrial = local_getNumericColumn(B, 'Trial');
        trialDistance = local_getNumericColumn(B, 'Vibrate_distance');
        monster = local_getNumericColumn(B, 'Monster');
        tactileOnset = local_getNumericColumn(B, 'Time_vibration');

        conditionTraces = cell(1, nConditions);

        for row = 1:height(B)
            conditionIndex = find(abs(analysisDistances - trialDistance(row)) < 1e-9, 1);

            % Only central trials with a tactile event are analysed.
            if monster(row) ~= 0 || isempty(conditionIndex)
                continue
            end

            trialID = behaviouralTrial(row);
            sampleMask = eyeTrial == trialID;

            retained = false;
            reason = '';
            pupilChange = NaN;
            baselineValidFraction = NaN;
            responseValidFraction = NaN;
            baselineSampleCount = 0;
            responseSampleCount = 0;

            if ~isfinite(tactileOnset(row))
                reason = 'Missing tactile onset';
            elseif ~any(sampleMask)
                reason = 'No matching eye samples for Trial ID';
            else
                relativeTime = eyeTimestamp(sampleMask) - tactileOnset(row);
                pupilValues = pupilAverage(sampleMask);

                baselineMask = relativeTime >= baselineWindow(1) & ...
                               relativeTime < baselineWindow(2);
                responseMask = relativeTime >= responseWindow(1) & ...
                               relativeTime <= responseWindow(2);

                baselineSampleCount = sum(baselineMask);
                responseSampleCount = sum(responseMask);

                if baselineSampleCount > 0
                    baselineValidFraction = sum(isfinite(pupilValues(baselineMask))) / ...
                                            baselineSampleCount;
                end
                if responseSampleCount > 0
                    responseValidFraction = sum(isfinite(pupilValues(responseMask))) / ...
                                            responseSampleCount;
                end

                if baselineSampleCount == 0 || responseSampleCount == 0
                    reason = 'Required time window absent';
                elseif baselineValidFraction < minimumValidFraction || ...
                       responseValidFraction < minimumValidFraction
                    reason = 'Insufficient valid pupil samples';
                else
                    baselineMean = local_meanOmitNaN(pupilValues(baselineMask));
                    responseMean = local_meanOmitNaN(pupilValues(responseMask));
                    pupilChange = responseMean - baselineMean;
                    retained = isfinite(pupilChange);

                    if retained
                        baselineCorrected = pupilValues - baselineMean;
                        trace = local_interpolateTrace(relativeTime, ...
                            baselineCorrected, commonTime);
                        conditionTraces{conditionIndex}(end+1,:) = trace; %#ok<AGROW>
                    else
                        reason = 'Invalid pupil-change estimate';
                    end
                end
            end

            trialRows(end+1,:) = {participant, phaseNames{p}, trialID, ...
                conditionLabels{conditionIndex}, trialDistance(row), pupilChange, ...
                baselineValidFraction, responseValidFraction, ...
                baselineSampleCount, responseSampleCount, retained, reason}; %#ok<AGROW>
        end

        for c = 1:nConditions
            if size(conditionTraces{c},1) >= minimumValidTrialsPerCondition
                meanTrace = local_meanOmitNaN(conditionTraces{c}, 1);
                participantTimeCourse(s,p,c,:) = reshape(meanTrace, 1, 1, 1, []);
            end
        end
    end
end

TrialResults = cell2table(trialRows, 'VariableNames', ...
    {'Subject','Phase','Trial','Condition','Distance','PupilChange_mm', ...
     'BaselineValidFraction','ResponseValidFraction','BaselineSamples', ...
     'ResponseSamples','Retained','ExclusionReason'});

TrialResults.Subject = string(TrialResults.Subject);
TrialResults.Phase = string(TrialResults.Phase);
TrialResults.Condition = string(TrialResults.Condition);
TrialResults.ExclusionReason = string(TrialResults.ExclusionReason);
writetable(TrialResults, fullfile(outputFolder, 'WP4_Pupil_Trial_Level_Results.csv'));

%% ================= PARTICIPANT-LEVEL SUMMARIES =======================

summaryRows = cell(0, 7);

for s = 1:nParticipants
    for p = 1:nPhases
        for c = 1:nConditions
            use = TrialResults.Subject == string(includedParticipantIDs{s}) & ...
                  TrialResults.Phase == string(phaseNames{p}) & ...
                  abs(TrialResults.Distance - analysisDistances(c)) < 1e-9 & ...
                  TrialResults.Retained;
            values = TrialResults.PupilChange_mm(use);
            validTrialCount = sum(isfinite(values));
            meetsMinimumTrials = validTrialCount >= minimumValidTrialsPerCondition;

            if meetsMinimumTrials
                participantConditionMean = local_meanOmitNaN(values);
            else
                participantConditionMean = NaN;
            end

            summaryRows(end+1,:) = {includedParticipantIDs{s}, phaseNames{p}, ...
                conditionLabels{c}, analysisDistances(c), ...
                participantConditionMean, validTrialCount, ...
                meetsMinimumTrials}; %#ok<AGROW>
        end
    end
end

ParticipantSummary = cell2table(summaryRows, 'VariableNames', ...
    {'Subject','Phase','Condition','Distance','MeanPupilChange_mm', ...
     'ValidTrials','MeetsMinimumTrials'});
ParticipantSummary.Subject = string(ParticipantSummary.Subject);
ParticipantSummary.Phase = string(ParticipantSummary.Phase);
ParticipantSummary.Condition = string(ParticipantSummary.Condition);
writetable(ParticipantSummary, ...
    fullfile(outputFolder, 'WP4_Pupil_Participant_Condition_Summary.csv'));

PupilCellQC = ParticipantSummary(:, ...
    {'Subject','Phase','Condition','Distance','ValidTrials', ...
     'MeetsMinimumTrials'});
writetable(PupilCellQC, ...
    fullfile(outputFolder, 'WP4_Pupil_Participant_Condition_QC.csv'));

% Matrices: participant x condition.
preValues = nan(nParticipants, nConditions);
postValues = nan(nParticipants, nConditions);

for s = 1:nParticipants
    for c = 1:nConditions
        preUse = ParticipantSummary.Subject == string(includedParticipantIDs{s}) & ...
                 ParticipantSummary.Phase == "Baseline" & ...
                 abs(ParticipantSummary.Distance - analysisDistances(c)) < 1e-9;
        postUse = ParticipantSummary.Subject == string(includedParticipantIDs{s}) & ...
                  ParticipantSummary.Phase == "Post-training" & ...
                  abs(ParticipantSummary.Distance - analysisDistances(c)) < 1e-9;

        preValues(s,c) = ParticipantSummary.MeanPupilChange_mm(preUse);
        postValues(s,c) = ParticipantSummary.MeanPupilChange_mm(postUse);
    end
end

%% ======================== EXACT INFERENCE ============================

conditionRows = cell(nConditions, 15);
for c = 1:nConditions
    pre = preValues(:,c);
    post = postValues(:,c);
    complete = isfinite(pre) & isfinite(post);
    pre = pre(complete);
    post = post(complete);
    change = post - pre;

    conditionRows(c,:) = {conditionLabels{c}, analysisDistances(c), numel(change), ...
        local_meanOmitNaN(pre), local_sampleSD(pre), ...
        local_meanOmitNaN(post), local_sampleSD(post), ...
        local_meanOmitNaN(change), local_sampleSD(change), ...
        local_exactSignFlip(pre), local_exactSignFlip(post), ...
        local_exactSignFlip(change), NaN, NaN, NaN};
end

ConditionResults = cell2table(conditionRows, 'VariableNames', ...
    {'Condition','Distance','N','BaselineMean_mm','BaselineSD_mm', ...
     'PostTrainingMean_mm','PostTrainingSD_mm', ...
     'PostTrainingMinusBaselineMean_mm', ...
     'PostTrainingMinusBaselineSD_mm', ...
     'BaselineVsZero_ExactP','PostTrainingVsZero_ExactP', ...
     'BaselinePostTraining_ExactP', ...
     'BaselineVsZero_BonferroniP', ...
     'PostTrainingVsZero_BonferroniP', ...
     'BaselinePostTraining_BonferroniP'});

ConditionResults.Condition = string(ConditionResults.Condition);
ConditionResults.BaselineVsZero_BonferroniP = ...
    min(1, ConditionResults.BaselineVsZero_ExactP * nConditions);
ConditionResults.PostTrainingVsZero_BonferroniP = ...
    min(1, ConditionResults.PostTrainingVsZero_ExactP * nConditions);
ConditionResults.BaselinePostTraining_BonferroniP = ...
    min(1, ConditionResults.BaselinePostTraining_ExactP * nConditions);
writetable(ConditionResults, ...
    fullfile(outputFolder, 'WP4_Pupil_Condition_Exact_Tests.csv'));

% Overall response averaged across the five visuo-tactile distances.
overallPre = local_meanOmitNaN(preValues(:,mainDistanceIndices), 2);
overallPost = local_meanOmitNaN(postValues(:,mainDistanceIndices), 2);
completeFiveDistanceProfiles = ...
    all(isfinite(preValues(:,mainDistanceIndices)),2) & ...
    all(isfinite(postValues(:,mainDistanceIndices)),2);
overallPre(~completeFiveDistanceProfiles) = NaN;
overallPost(~completeFiveDistanceProfiles) = NaN;
overallChange = overallPost - overallPre;
overallP = local_exactSignFlip(overallChange);

% Distance-response slopes across the five visuo-tactile distances.
preSlopePupil = nan(nParticipants,1);
postSlopePupil = nan(nParticipants,1);
for s = 1:nParticipants
    yPre = preValues(s,mainDistanceIndices);
    yPost = postValues(s,mainDistanceIndices);
    if all(isfinite(yPre))
        fitPre = polyfit(ppsDistances, yPre, 1);
        preSlopePupil(s) = fitPre(1);
    end
    if all(isfinite(yPost))
        fitPost = polyfit(ppsDistances, yPost, 1);
        postSlopePupil(s) = fitPost(1);
    end
end
slopeChange = postSlopePupil - preSlopePupil;
slopeP = local_exactSignFlip(slopeChange);

OverallResults = table( ...
    ["Mean response across five PPS distances"; "Distance-response slope"], ...
    [sum(isfinite(overallChange)); sum(isfinite(slopeChange))], ...
    [local_meanOmitNaN(overallPre); local_meanOmitNaN(preSlopePupil)], ...
    [local_sampleSD(overallPre); local_sampleSD(preSlopePupil)], ...
    [local_meanOmitNaN(overallPost); local_meanOmitNaN(postSlopePupil)], ...
    [local_sampleSD(overallPost); local_sampleSD(postSlopePupil)], ...
    [local_meanOmitNaN(overallChange); local_meanOmitNaN(slopeChange)], ...
    [local_sampleSD(overallChange); local_sampleSD(slopeChange)], ...
    [overallP; slopeP], ...
    'VariableNames', {'Outcome','N','BaselineMean','BaselineSD', ...
                      'PostTrainingMean','PostTrainingSD', ...
                      'PostTrainingMinusBaselineMean', ...
                      'PostTrainingMinusBaselineSD','ExactP'});
writetable(OverallResults, ...
    fullfile(outputFolder, 'WP4_Pupil_Overall_Exact_Tests.csv'));

%% ============================== FIGURES ==============================

preColor = [0.90 0.12 0.12];
postColor = [0.12 0.28 0.90];
lightPre = [1.00 0.72 0.72];
lightPost = [0.72 0.78 1.00];

% Figure 1: WP3-style distance profile, shown for Baseline and Post-training.
figure(1); clf; hold on
x = ppsDistances;
for s = 1:nParticipants
    plot(x, preValues(s,mainDistanceIndices), '-', 'Color', lightPre, ...
        'LineWidth', 0.8, 'HandleVisibility','off');
    plot(x, postValues(s,mainDistanceIndices), '-', 'Color', lightPost, ...
        'LineWidth', 0.8, 'HandleVisibility','off');
end

preMean = local_meanOmitNaN(preValues(:,mainDistanceIndices), 1);
postMean = local_meanOmitNaN(postValues(:,mainDistanceIndices), 1);
preSEM = local_sem(preValues(:,mainDistanceIndices), 1);
postSEM = local_sem(postValues(:,mainDistanceIndices), 1);

errorbar(x, preMean, preSEM, '-o', 'Color', preColor, ...
    'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',2.2, ...
    'DisplayName','BASELINE');
errorbar(x, postMean, postSEM, '-o', 'Color', postColor, ...
    'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',2.2, ...
    'DisplayName','POST-TRAINING');
yline(0, 'k--', 'HandleVisibility','off');
xlabel('Distance (m)');
ylabel('\Delta pupil diameter (mm)');
xticks(ppsDistances);
xticklabels(compose('%.2f', ppsDistances));
legend('Location','best');
box off; set(gca,'TickDir','out','LineWidth',1.5,'FontSize',13);
set(gcf,'Color','w','Position',[100 100 650 520]);
exportgraphics(gcf, fullfile(outputFolder, ...
    'Figure_1_Pupil_Change_by_Distance_Baseline_PostTraining.png'), 'Resolution',300);
savefig(gcf, fullfile(outputFolder, ...
    'Figure_1_Pupil_Change_by_Distance_Baseline_PostTraining.fig'));

% Figure 2: paired overall response across five visuo-tactile distances.
figure(2); clf; hold on
for s = 1:nParticipants
    plot([1 2], [overallPre(s) overallPost(s)], '-', ...
        'Color',[0.72 0.72 0.72], 'LineWidth',1.2);
end
scatter(ones(nParticipants,1), overallPre, 55, 'o', ...
    'MarkerEdgeColor',preColor, 'MarkerFaceColor','w', 'LineWidth',1.5);
scatter(2*ones(nParticipants,1), overallPost, 55, 'o', ...
    'MarkerEdgeColor',postColor, 'MarkerFaceColor','w', 'LineWidth',1.5);
errorbar(1, local_meanOmitNaN(overallPre), local_sem(overallPre,1), 'o', ...
    'Color',preColor,'MarkerFaceColor',preColor,'MarkerSize',9,'LineWidth',2);
errorbar(2, local_meanOmitNaN(overallPost), local_sem(overallPost,1), 'o', ...
    'Color',postColor,'MarkerFaceColor',postColor,'MarkerSize',9,'LineWidth',2);
yline(0, 'k--', 'HandleVisibility','off');
xlim([0.6 2.4]); xticks([1 2]);
xticklabels({'BASELINE','POST-TRAINING'});
ylabel('Mean \Delta pupil diameter (mm)');
title(sprintf('Five-distance mean: exact p = %.4f', overallP));
box off; set(gca,'TickDir','out','LineWidth',1.5,'FontSize',13);
set(gcf,'Color','w','Position',[120 120 520 520]);
exportgraphics(gcf, fullfile(outputFolder, ...
    'Figure_2_Overall_Pupil_Response_Baseline_PostTraining.png'), 'Resolution',300);
savefig(gcf, fullfile(outputFolder, ...
    'Figure_2_Overall_Pupil_Response_Baseline_PostTraining.fig'));

% Figure 3: event-related pupil time courses for all tactile conditions.
figure(3); clf
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for c = 1:nConditions
    nexttile; hold on
    preCourse = squeeze(participantTimeCourse(:,1,c,:));
    postCourse = squeeze(participantTimeCourse(:,2,c,:));
    local_plotMeanSEM(commonTime, preCourse, preColor);
    local_plotMeanSEM(commonTime, postCourse, postColor);
    xline(0,'k--'); yline(0,'k:');
    xlim([baselineWindow(1) responseWindow(2)]);
    title(conditionLabels{c});
    xlabel('Time from tactile onset (s)');
    ylabel('\Delta pupil diameter (mm)');
    box off; set(gca,'TickDir','out','LineWidth',1.2,'FontSize',10);
    if c == 1
        plot(nan,nan,'-','Color',preColor,'LineWidth',2, ...
            'DisplayName','BASELINE');
        plot(nan,nan,'-','Color',postColor,'LineWidth',2, ...
            'DisplayName','POST-TRAINING');
        legend('Location','best');
    end
end
set(gcf,'Color','w','Position',[80 80 1100 680]);
exportgraphics(gcf, fullfile(outputFolder, ...
    'Figure_3_Event_Related_Pupil_Time_Courses.png'), 'Resolution',300);
savefig(gcf, fullfile(outputFolder, ...
    'Figure_3_Event_Related_Pupil_Time_Courses.fig'));

%% ============================ REPORTING ==============================

fprintf('\n============================================================\n');
fprintf('WP4 PUPILLOMETRY PRELIMINARY RESULTS\n');
fprintf('============================================================\n');
fprintf('Participants included: %d\n', nParticipants);
fprintf('Baseline window: %.1f to %.1f s relative to tactile onset\n', ...
    baselineWindow(1), baselineWindow(2));
fprintf('Response window: %.1f to %.1f s relative to tactile onset\n', ...
    responseWindow(1), responseWindow(2));
fprintf('Minimum valid pupil fraction per window: %.0f%%\n\n', ...
    100*minimumValidFraction);
fprintf('Minimum valid trials per participant/phase/condition: %d of 16\n\n', ...
    minimumValidTrialsPerCondition);

fprintf('Overall response across five visuo-tactile distances:\n');
fprintf('  Baseline:      M = %.4f, SD = %.4f mm\n', ...
    local_meanOmitNaN(overallPre), local_sampleSD(overallPre));
fprintf('  Post-training: M = %.4f, SD = %.4f mm\n', ...
    local_meanOmitNaN(overallPost), local_sampleSD(overallPost));
fprintf('  Post-training - baseline: M = %.4f, SD = %.4f mm\n', ...
    local_meanOmitNaN(overallChange), local_sampleSD(overallChange));
fprintf('  Exact two-sided sign-flip p = %.5f\n\n', overallP);

fprintf('Distance-response slope:\n');
fprintf('  Baseline:      M = %.4f mm/m\n', local_meanOmitNaN(preSlopePupil));
fprintf('  Post-training: M = %.4f mm/m\n', local_meanOmitNaN(postSlopePupil));
fprintf('  Post-training - baseline: M = %.4f mm/m\n', ...
    local_meanOmitNaN(slopeChange));
fprintf('  Exact two-sided sign-flip p = %.5f\n\n', slopeP);

fprintf('Condition-specific baseline/Post-training comparisons:\n');
disp(ConditionResults(:, {'Condition','N','BaselineMean_mm', ...
    'PostTrainingMean_mm','PostTrainingMinusBaselineMean_mm', ...
    'BaselinePostTraining_ExactP', ...
    'BaselinePostTraining_BonferroniP'}));

fprintf('\nParticipant/phase/condition cells below the valid-trial threshold:\n');
disp(PupilCellQC(~PupilCellQC.MeetsMinimumTrials,:));

fprintf(['\nAll pupil findings are exploratory and preliminary. ', ...
         'The visual-only (-1) condition was not included because ', ...
         'there is no tactile onset for tactile-event alignment.\n']);
fprintf('All outputs were saved to:\n%s\n', outputFolder);
fprintf('============================================================\n');

save(fullfile(outputFolder, 'WP4_Pupillometry_Preliminary_Workspace.mat'));

%% =========================== LOCAL FUNCTIONS =========================

function filePath = local_resolveFile(folderPath, stem, fileRole)
    if exist(folderPath, 'dir') ~= 7
        error('The configured %s folder does not exist:\n%s', ...
            fileRole, folderPath);
    end

    if strcmpi(fileRole, 'eye')
        exactCandidates = {[stem '.csv'], ['eye_' stem '.csv']};
        wildcardPatterns = {[stem '*.csv'], ['eye_' stem '*.csv']};
    else
        exactCandidates = {[stem '(1).csv'], [stem '_behaviour.csv'], ...
                           [stem '_behavior.csv'], [stem '.csv']};
        wildcardPatterns = {[stem '*.csv'], ...
                            [stem '*behaviour*.csv'], ...
                            [stem '*behavior*.csv']};
    end

    % Try the expected exact names first, then any same-stem CSV. Every
    % candidate is classified by its header, so behavioural and eye files
    % can safely coexist in the same folder even when Windows adds suffixes.
    candidatePaths = cell(0,1);

    for i = 1:numel(exactCandidates)
        candidatePaths{end+1,1} = fullfile( ...
            folderPath, exactCandidates{i}); %#ok<AGROW>
    end

    for patternNumber = 1:numel(wildcardPatterns)
        listing = dir(fullfile(folderPath, ...
            wildcardPatterns{patternNumber}));
        for fileNumber = 1:numel(listing)
            if ~listing(fileNumber).isdir
                candidatePaths{end+1,1} = fullfile( ...
                    listing(fileNumber).folder, ...
                    listing(fileNumber).name); %#ok<AGROW>
            end
        end
    end

    if ~isempty(candidatePaths)
        candidatePaths = unique(candidatePaths, 'stable');
    end

    matchingPaths = cell(0,1);
    for i = 1:numel(candidatePaths)
        candidate = candidatePaths{i};
        if exist(candidate, 'file') == 2 && ...
                local_fileMatchesRole(candidate, fileRole)
            matchingPaths{end+1,1} = candidate; %#ok<AGROW>
        end
    end

    if numel(matchingPaths) == 1
        filePath = matchingPaths{1};
        fprintf('Resolved %s file for %s:\n  %s\n', ...
            fileRole, stem, filePath);
        return
    end

    if numel(matchingPaths) > 1
        error(['More than one %s file matched %s. Move duplicates out ' ...
               'of the folder or rename the intended file:\n%s'], ...
            fileRole, stem, strjoin(matchingPaths, '\n'));
    end

    availableFiles = dir(fullfile(folderPath, '*.csv'));
    availableNames = string({availableFiles.name});
    relatedNames = availableNames(contains( ...
        lower(availableNames), lower(string(stem))));

    if isempty(relatedNames)
        relatedText = '(no CSV filenames containing this participant/phase stem)';
    else
        relatedText = strjoin(cellstr(relatedNames(:)), ', ');
    end

    error(['Could not find a %s CSV for %s in:\n%s\n' ...
           'Files with a related name: %s\n' ...
           'Check the folder setting and confirm that the file header ' ...
           'contains the required %s columns.'], ...
        fileRole, stem, folderPath, relatedText, fileRole);
end

function matches = local_fileMatchesRole(filePath, fileRole)
    fileID = fopen(filePath, 'r');
    if fileID < 0
        matches = false;
        return
    end
    cleanUp = onCleanup(@() fclose(fileID)); %#ok<NASGU>
    header = lower(fgetl(fileID));

    if strcmpi(fileRole, 'eye')
        matches = contains(header, 'timestamp') && ...
                  contains(header, 'pupil_diameter');
    else
        matches = contains(header, 'time_vibration') && ...
                  contains(header, 'vibrate_distance') && ...
                  contains(header, 'monster');
    end
end

function T = local_readAndCleanNames(filePath)
    T = readtable(filePath, 'VariableNamingRule','preserve');
    names = T.Properties.VariableNames;
    for i = 1:numel(names)
        names{i} = strtrim(names{i});
    end
    T.Properties.VariableNames = names;
end

function values = local_getNumericColumn(T, requestedName)
    names = T.Properties.VariableNames;
    normalisedNames = cellfun(@local_normaliseName, names, ...
        'UniformOutput',false);
    requested = local_normaliseName(requestedName);
    index = find(strcmp(normalisedNames, requested), 1);

    if isempty(index)
        error('Required column "%s" was not found.', requestedName);
    end

    raw = T.(names{index});
    if isnumeric(raw) || islogical(raw)
        values = double(raw);
    else
        values = str2double(string(raw));
    end
    values = values(:);
end

function name = local_normaliseName(name)
    name = lower(regexprep(strtrim(char(name)), '[^a-zA-Z0-9]', ''));
end

function average = local_rowMeanTwoEyes(left, right)
    average = nan(size(left));
    both = isfinite(left) & isfinite(right);
    leftOnly = isfinite(left) & ~isfinite(right);
    rightOnly = ~isfinite(left) & isfinite(right);
    average(both) = (left(both) + right(both)) / 2;
    average(leftOnly) = left(leftOnly);
    average(rightOnly) = right(rightOnly);
end

function trace = local_interpolateTrace(time, values, commonTime)
    valid = isfinite(time) & isfinite(values);
    time = time(valid);
    values = values(valid);

    if numel(time) < 2
        trace = nan(size(commonTime));
        return
    end

    [time, uniqueIndex] = unique(time, 'stable');
    values = values(uniqueIndex);
    if numel(time) < 2
        trace = nan(size(commonTime));
        return
    end

    trace = interp1(time, values, commonTime, 'linear', NaN);
end

function value = local_meanOmitNaN(x, dim)
    if nargin < 2
        x = x(:);
        valid = isfinite(x);
        if ~any(valid)
            value = NaN;
        else
            value = mean(x(valid));
        end
        return
    end

    valid = isfinite(x);
    x(~valid) = 0;
    counts = sum(valid, dim);
    value = sum(x, dim) ./ counts;
    value(counts == 0) = NaN;
end

function value = local_sampleSD(x)
    x = x(isfinite(x));
    if numel(x) < 2
        value = NaN;
    else
        value = sqrt(sum((x - mean(x)).^2) / (numel(x)-1));
    end
end

function sem = local_sem(x, dim)
    if nargin < 2, dim = 1; end
    valid = isfinite(x);
    n = sum(valid, dim);
    meanValue = local_meanOmitNaN(x, dim);

    difference = x - meanValue;
    difference(~valid) = 0;
    variance = sum(difference.^2, dim) ./ max(n-1, 1);
    sd = sqrt(variance);
    sem = sd ./ sqrt(n);
    sem(n < 2) = NaN;
end

function p = local_exactSignFlip(values)
    values = values(isfinite(values));
    values = values(:);
    n = numel(values);

    if n == 0
        p = NaN;
        return
    end

    observed = abs(mean(values));
    nPermutations = 2^n;
    permutationStatistics = nan(nPermutations,1);

    for permutation = 0:(nPermutations-1)
        bits = bitget(uint64(permutation), 1:n);
        signs = 2*double(bits(:)) - 1;
        permutationStatistics(permutation+1) = abs(mean(values .* signs));
    end

    tolerance = 1e-12;
    p = mean(permutationStatistics >= observed - tolerance);
end

function local_plotMeanSEM(time, participantCourses, color)
    if isvector(participantCourses)
        participantCourses = participantCourses(:)';
    end
    meanCourse = local_meanOmitNaN(participantCourses, 1);
    semCourse = local_sem(participantCourses, 1);
    valid = isfinite(meanCourse) & isfinite(semCourse);

    if any(valid)
        x = time(valid);
        upper = meanCourse(valid) + semCourse(valid);
        lower = meanCourse(valid) - semCourse(valid);
        fill([x fliplr(x)], [upper fliplr(lower)], color, ...
            'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off');
        plot(x, meanCourse(valid), 'Color',color,'LineWidth',2, ...
            'HandleVisibility','off');
    end
end

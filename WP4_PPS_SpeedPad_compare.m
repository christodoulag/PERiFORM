%% WP4 INTEGRATED PPS AND SPEEDPAD PRELIMINARY ANALYSIS
% Standalone analysis of the principal WP4 outcomes.
%
% It performs:
%   1. Baseline PPS quality control using the established WP3 criteria.
%   2. Participant-level baseline and Post-training PPS-border estimation.
%   3. Exact paired sign-flip analysis of PPS Post-training - baseline
%      change.
%   4. SpeedPad Time-by-Condition interaction contrast:
%        (PULL TRAINING Second - First) - (CONTROL Second - First).
%   5. Exact Spearman association between PPS change and SpeedPad benefit.
%
% No Statistics and Machine Learning Toolbox functions are required.
% The sigmoid fitting step uses FIT/FITTYPE, as in the existing WP3 code.

clear; close all; clc;
dbstop if error
rng(20260804,'twister');

%% ======================== USER SETTINGS ==============================

% Folder containing the behavioural PPS CSV files.
ppsDataFolder = '';

% SpeedPad workbook and worksheet.
speedPadFile = fullfile(ppsDataFolder,'speedpabdata.xlsx');
speedPadSheet = 'Sheet1';

% Analysis outputs are placed in a separate subfolder.
outputFolder = fullfile(ppsDataFolder,'WP4_Integrated_Outputs');

% Participant identifiers corresponding to the PPS filenames.
participantNames = {''};

distances = [0.25; 0.75; 1.25; 1.75; 2.25];

% Established PPS cleaning and baseline-QC parameters.
minimumRT = 0.1;
maximumRT = 1.0;
outlierCutoff = 3;
minimumBaselineR2 = 0.10;

% This threshold is reported as a fit-quality flag only. Post-training data
% are not excluded on this basis because doing so after Pull Training could
% bias the estimated change. The flag supports transparent interpretation.
fitQualityFlagThreshold = 0.80;

makeFigures = true;

if ~isfolder(ppsDataFolder)
    error('The PPS data folder was not found:\n%s',ppsDataFolder);
end
if ~isfile(speedPadFile)
    error('The SpeedPad workbook was not found:\n%s',speedPadFile);
end
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

%% ======================== SIGMOID MODEL ==============================

eqsig = fittype( ...
    '1/(1+exp(-(x-x0)/b))', ...
    'independent','x', ...
    'coefficients',{'b','x0'});

fitOptions = fitoptions(eqsig);
fitOptions.Lower = [eps 0.25];
fitOptions.Upper = [Inf 2.25];
fitOptions.StartPoint = [0.5 1.25];

%% ======================== PPS ANALYSIS ===============================

nCandidateParticipants = numel(participantNames);
nDistances = numel(distances);

numericParticipantID = nan(nCandidateParticipants,1);
preMeanRT = nan(nCandidateParticipants,nDistances);
postMeanRT = nan(nCandidateParticipants,nDistances);
preValidTrials = nan(nCandidateParticipants,nDistances);
postValidTrials = nan(nCandidateParticipants,nDistances);

preLinearSlope = nan(nCandidateParticipants,1);
preLinearR2 = nan(nCandidateParticipants,1);
prePPSBorder = nan(nCandidateParticipants,1);
postPPSBorder = nan(nCandidateParticipants,1);
preSigmoidR2 = nan(nCandidateParticipants,1);
postSigmoidR2 = nan(nCandidateParticipants,1);

for s = 1:nCandidateParticipants
    numericParticipantID(s) = local_numericParticipantID(participantNames{s});

    preFile = local_resolvePPSFile( ...
        ppsDataFolder,participantNames{s},'pre');
    postFile = local_resolvePPSFile( ...
        ppsDataFolder,participantNames{s},'post');

    [preMeanRT(s,:),preValidTrials(s,:)] = ...
        local_conditionMeans(preFile,distances,minimumRT, ...
        maximumRT,outlierCutoff);

    [postMeanRT(s,:),postValidTrials(s,:)] = ...
        local_conditionMeans(postFile,distances,minimumRT, ...
        maximumRT,outlierCutoff);

    % Baseline QC uses the same linear slope and R2 logic as WP3.
    if all(isfinite(preMeanRT(s,:))) && ...
            max(preMeanRT(s,:)) > min(preMeanRT(s,:))
        linearCoefficients = polyfit(distances,preMeanRT(s,:)',1);
        fittedLinear = polyval(linearCoefficients,distances);
        residualSS = sum((preMeanRT(s,:)'-fittedLinear).^2);
        totalSS = sum((preMeanRT(s,:)'-mean(preMeanRT(s,:))).^2);
        preLinearSlope(s) = linearCoefficients(1);
        if totalSS > 0
            preLinearR2(s) = 1-residualSS/totalSS;
        end
    end

    [prePPSBorder(s),preSigmoidR2(s)] = local_fitPPSBorder( ...
        distances,preMeanRT(s,:)',eqsig,fitOptions);

    [postPPSBorder(s),postSigmoidR2(s)] = local_fitPPSBorder( ...
        distances,postMeanRT(s,:)',eqsig,fitOptions);
end

lowOrInvalidPreR2 = ...
    ~isfinite(preLinearR2) | preLinearR2 < minimumBaselineR2;
nonPositivePreSlope = ...
    ~isfinite(preLinearSlope) | preLinearSlope <= 0;
excludedBaselineQC = lowOrInvalidPreR2 | nonPositivePreSlope;

invalidPPSFit = ~isfinite(prePPSBorder) | ~isfinite(postPPSBorder);
includedPPS = ~excludedBaselineQC & ~invalidPPSFit;

preFitFlag = preSigmoidR2 < fitQualityFlagThreshold | ...
             ~isfinite(preSigmoidR2);
postFitFlag = postSigmoidR2 < fitQualityFlagThreshold | ...
              ~isfinite(postSigmoidR2);

exclusionReason = repmat("Included",nCandidateParticipants,1);
for s = 1:nCandidateParticipants
    if invalidPPSFit(s)
        exclusionReason(s) = "Invalid baseline or Post-training sigmoid fit";
    elseif lowOrInvalidPreR2(s) && nonPositivePreSlope(s)
        exclusionReason(s) = "Baseline R2 < 0.10 and non-positive slope";
    elseif lowOrInvalidPreR2(s)
        exclusionReason(s) = "Baseline R2 < 0.10";
    elseif nonPositivePreSlope(s)
        exclusionReason(s) = "Non-positive baseline slope";
    end
end

PPSSummaryAll = table( ...
    numericParticipantID, ...
    string(participantNames(:)), ...
    preLinearSlope,preLinearR2, ...
    prePPSBorder,preSigmoidR2,preFitFlag, ...
    postPPSBorder,postSigmoidR2,postFitFlag, ...
    excludedBaselineQC,invalidPPSFit,includedPPS,exclusionReason, ...
    'VariableNames',{ ...
    'ParticipantID','ParticipantCode','BaselineLinearSlope', ...
    'BaselineLinearR2','BaselinePPSBorder_m','BaselineSigmoidR2', ...
    'BaselineFitR2Below_0_80','PostTrainingPPSBorder_m', ...
    'PostTrainingSigmoidR2','PostTrainingFitR2Below_0_80', ...
    'ExcludedBaselineQC','InvalidPPSFit','IncludedIntegrated', ...
    'ExclusionReason'});

writetable(PPSSummaryAll,fullfile(outputFolder, ...
    'WP4_Integrated_PPS_QC_and_Fit_Summary.csv'));

fprintf('\nPPS baseline QC and participant-level fitting completed.\n');
fprintf('Candidate participants: %d\n',nCandidateParticipants);
fprintf('Included participants:  %d\n',sum(includedPPS));
if any(~includedPPS)
    fprintf('Excluded participant(s):\n');
    disp(PPSSummaryAll(~includedPPS, ...
        {'ParticipantCode','BaselineLinearSlope','BaselineLinearR2', ...
         'ExclusionReason'}));
end

%% ======================== SPEEDPAD IMPORT ============================

rawSpeedPad = readcell(speedPadFile,'Sheet',speedPadSheet);

if size(rawSpeedPad,1) < 2 || size(rawSpeedPad,2) < 14
    error(['The SpeedPad workbook does not match the expected layout. ' ...
        'Expected data in columns A:N.']);
end

speedPadCandidateIDs = local_extractCounts( ...
    rawSpeedPad(2:end,1),false,'SN');
speedPadRowMask = isfinite(speedPadCandidateIDs);

if ~any(speedPadRowMask)
    error('No SpeedPad participant rows were detected in column A.');
end

speedPadIDs = speedPadCandidateIDs(speedPadRowMask);
speedPadRows = rawSpeedPad(1+find(speedPadRowMask),:);

pullTrainingLabels1 = strtrim(string(speedPadRows(:,2)));
pullTrainingLabels2 = strtrim(string(speedPadRows(:,5)));
controlLabels1 = strtrim(string(speedPadRows(:,9)));
controlLabels2 = strtrim(string(speedPadRows(:,12)));

% The existing workbook may contain the historical label "Intervention".
% It is accepted as an input alias, but all new outputs use PULL TRAINING.
validPullTrainingLabels1 = strcmpi(pullTrainingLabels1,'Intervention') | ...
    strcmpi(pullTrainingLabels1,'Pull Training') | ...
    strcmpi(pullTrainingLabels1,'PullTraining');
validPullTrainingLabels2 = strcmpi(pullTrainingLabels2,'Intervention') | ...
    strcmpi(pullTrainingLabels2,'Pull Training') | ...
    strcmpi(pullTrainingLabels2,'PullTraining');

if ~all(validPullTrainingLabels1) || ~all(validPullTrainingLabels2)
    error(['Pull Training labels were not found in columns B and E. ' ...
        'Accepted input labels are Pull Training, PullTraining, ' ...
        'or the historical label Intervention.']);
end
if ~all(strcmpi(controlLabels1,'Control')) || ...
        ~all(strcmpi(controlLabels2,'Control'))
    error('Control labels were not found in columns I and L.');
end

pullTrainingFirstAssessment = local_extractCounts( ...
    speedPadRows(:,3),true,'Pull Training First Assessment Hits');
pullTrainingSecondAssessment = local_extractCounts( ...
    speedPadRows(:,6),true,'Pull Training Second Assessment Hits');
controlFirstAssessment = local_extractCounts( ...
    speedPadRows(:,10),true,'Control First Assessment Hits');
controlSecondAssessment = local_extractCounts( ...
    speedPadRows(:,13),true,'Control Second Assessment Hits');

if numel(unique(speedPadIDs)) ~= numel(speedPadIDs)
    error('SpeedPad participant identifiers must be unique.');
end

SpeedPadAll = table(speedPadIDs,pullTrainingFirstAssessment, ...
    pullTrainingSecondAssessment,controlFirstAssessment, ...
    controlSecondAssessment, ...
    'VariableNames',{'ParticipantID','PullTraining_FirstAssessment', ...
    'PullTraining_SecondAssessment','Control_FirstAssessment', ...
    'Control_SecondAssessment'});

SpeedPadAll.PullTrainingChange = ...
    SpeedPadAll.PullTraining_SecondAssessment- ...
    SpeedPadAll.PullTraining_FirstAssessment;
SpeedPadAll.ControlChange = ...
    SpeedPadAll.Control_SecondAssessment-SpeedPadAll.Control_FirstAssessment;
SpeedPadAll.SpeedPadBenefit = ...
    SpeedPadAll.PullTrainingChange-SpeedPadAll.ControlChange;

%% ======================== MATCH PARTICIPANTS =========================

PPSIncluded = PPSSummaryAll(includedPPS,:);

[matched,locationInSpeedPad] = ismember( ...
    PPSIncluded.ParticipantID,SpeedPadAll.ParticipantID);

if any(~matched)
    error('Missing SpeedPad data for participant ID(s): %s', ...
        mat2str(PPSIncluded.ParticipantID(~matched)'));
end

IntegratedData = PPSIncluded(:,{ ...
    'ParticipantID','ParticipantCode','BaselinePPSBorder_m', ...
    'BaselineSigmoidR2','BaselineFitR2Below_0_80', ...
    'PostTrainingPPSBorder_m','PostTrainingSigmoidR2', ...
    'PostTrainingFitR2Below_0_80'});

matchedSpeedPad = SpeedPadAll(locationInSpeedPad,:);
IntegratedData.PullTraining_FirstAssessment = ...
    matchedSpeedPad.PullTraining_FirstAssessment;
IntegratedData.PullTraining_SecondAssessment = ...
    matchedSpeedPad.PullTraining_SecondAssessment;
IntegratedData.PullTrainingChange = matchedSpeedPad.PullTrainingChange;
IntegratedData.Control_FirstAssessment = matchedSpeedPad.Control_FirstAssessment;
IntegratedData.Control_SecondAssessment = matchedSpeedPad.Control_SecondAssessment;
IntegratedData.ControlChange = matchedSpeedPad.ControlChange;
IntegratedData.PPSChange_m = ...
    IntegratedData.PostTrainingPPSBorder_m- ...
    IntegratedData.BaselinePPSBorder_m;
IntegratedData.SpeedPadBenefit_hits = matchedSpeedPad.SpeedPadBenefit;

nParticipants = height(IntegratedData);
if nParticipants < 2
    error('At least two matched participants are required.');
end

writetable(IntegratedData,fullfile(outputFolder, ...
    'WP4_Integrated_Participant_Data.csv'));

%% ======================== EXACT MAIN TESTS ===========================

[ppsMeanChange,ppsPtwo,ppsPgreater,ppsPless, ...
    ppsPermutations,ppsTestMethod] = ...
    local_signFlip(IntegratedData.PPSChange_m);

[speedMeanBenefit,speedPtwo,speedPgreater,speedPless, ...
    speedPermutations,speedTestMethod] = ...
    local_signFlip(IntegratedData.SpeedPadBenefit_hits);

MainTests = table( ...
    ["PPS Post-training - baseline"; ...
     "SpeedPad Pull Training change - Control change"], ...
    ["m";"hits"], ...
    [nParticipants;nParticipants], ...
    [mean(IntegratedData.BaselinePPSBorder_m); ...
     mean(IntegratedData.ControlChange)], ...
    [std(IntegratedData.BaselinePPSBorder_m); ...
     std(IntegratedData.ControlChange)], ...
    [mean(IntegratedData.PostTrainingPPSBorder_m); ...
     mean(IntegratedData.PullTrainingChange)], ...
    [std(IntegratedData.PostTrainingPPSBorder_m); ...
     std(IntegratedData.PullTrainingChange)], ...
    [ppsMeanChange;speedMeanBenefit], ...
    [std(IntegratedData.PPSChange_m); ...
     std(IntegratedData.SpeedPadBenefit_hits)], ...
    [ppsPtwo;speedPtwo], ...
    [ppsPgreater;speedPgreater], ...
    [ppsPless;speedPless], ...
    [ppsPermutations;speedPermutations], ...
    [ppsTestMethod;speedTestMethod], ...
    'VariableNames',{'Outcome','Unit','N','MeasureA_Mean','MeasureA_SD', ...
    'MeasureB_Mean','MeasureB_SD','ContrastMean','ContrastSD', ...
    'ExactP_TwoSided','ExactP_GreaterThanZero','ExactP_LessThanZero', ...
    'NumberOfPermutations','PermutationMethod'});

writetable(MainTests,fullfile(outputFolder, ...
    'WP4_Integrated_Main_Exact_Tests.csv'));

%% ======================== PPS-SPEEDPAD LINK ==========================

[spearmanRho,spearmanP,nSpearmanPermutations,spearmanMethod] = ...
    local_exactSpearman( ...
    IntegratedData.PPSChange_m, ...
    IntegratedData.SpeedPadBenefit_hits);

AssociationResults = table( ...
    nParticipants,spearmanRho,spearmanP,nSpearmanPermutations, ...
    spearmanMethod, ...
    'VariableNames',{'N','SpearmanRho','ExactP_TwoSided', ...
    'NumberOfPermutations','PermutationMethod'});

writetable(AssociationResults,fullfile(outputFolder, ...
    'WP4_PPS_Change_SpeedPad_Benefit_Association.csv'));

%% ======================== FIGURES ====================================

if makeFigures
    red = [0.90 0.12 0.12];
    blue = [0.12 0.28 0.90];
    grey = [0.72 0.72 0.72];

    figIntegrated = figure(1); clf(figIntegrated);
    set(figIntegrated,'Color','w','Position',[70 100 1400 480]);
    layout = tiledlayout(figIntegrated,1,3, ...
        'TileSpacing','compact','Padding','compact');

    % A. Participant-level PPS borders.
    axPPS = nexttile(layout,1); hold(axPPS,'on');
    for s = 1:nParticipants
        plot(axPPS,[1 2], ...
            [IntegratedData.BaselinePPSBorder_m(s), ...
             IntegratedData.PostTrainingPPSBorder_m(s)], ...
            '-','Color',grey,'LineWidth',1.1,'HandleVisibility','off');
    end
    scatter(axPPS,ones(nParticipants,1), ...
        IntegratedData.BaselinePPSBorder_m,55,'o', ...
        'MarkerEdgeColor',red,'MarkerFaceColor','w','LineWidth',1.4);
    scatter(axPPS,2*ones(nParticipants,1), ...
        IntegratedData.PostTrainingPPSBorder_m,55,'o', ...
        'MarkerEdgeColor',blue,'MarkerFaceColor','w','LineWidth',1.4);
    local_groupMeanSEM(axPPS,1,IntegratedData.BaselinePPSBorder_m,red);
    local_groupMeanSEM(axPPS,2,IntegratedData.PostTrainingPPSBorder_m,blue);
    xlim(axPPS,[0.6 2.4]);
    xticks(axPPS,[1 2]);
    xticklabels(axPPS,{'BASELINE','POST-TRAINING'});
    ylabel(axPPS,'PPS border (m)');
    title(axPPS,sprintf('PPS change: exact p = %.3f',ppsPtwo));
    local_styleAxes(axPPS);

    % B. Control and Pull Training changes in SpeedPad hits.
    axSpeed = nexttile(layout,2); hold(axSpeed,'on');
    for s = 1:nParticipants
        plot(axSpeed,[1 2], ...
            [IntegratedData.ControlChange(s), ...
             IntegratedData.PullTrainingChange(s)], ...
            '-','Color',grey,'LineWidth',1.1,'HandleVisibility','off');
    end
    scatter(axSpeed,ones(nParticipants,1), ...
        IntegratedData.ControlChange,55,'o', ...
        'MarkerEdgeColor',blue,'MarkerFaceColor','w','LineWidth',1.4);
    scatter(axSpeed,2*ones(nParticipants,1), ...
        IntegratedData.PullTrainingChange,55,'o', ...
        'MarkerEdgeColor',red,'MarkerFaceColor','w','LineWidth',1.4);
    local_groupMeanSEM(axSpeed,1,IntegratedData.ControlChange,blue);
    local_groupMeanSEM(axSpeed,2,IntegratedData.PullTrainingChange,red);
    yline(axSpeed,0,'k--','HandleVisibility','off');
    xlim(axSpeed,[0.6 2.4]);
    xticks(axSpeed,[1 2]);
    xticklabels(axSpeed,{'CONTROL','PULL TRAINING'});
    ylabel(axSpeed,'Change in hits');
    title(axSpeed,sprintf('SpeedPad benefit: exact p = %.3f',speedPtwo));
    local_styleAxes(axSpeed);

    % C. PPS change and Pull Training-specific SpeedPad benefit.
    axLink = nexttile(layout,3); hold(axLink,'on');
    scatter(axLink,IntegratedData.PPSChange_m, ...
        IntegratedData.SpeedPadBenefit_hits,70,'o', ...
        'MarkerEdgeColor',blue,'MarkerFaceColor','w','LineWidth',1.5);

    xValues = IntegratedData.PPSChange_m;
    yValues = IntegratedData.SpeedPadBenefit_hits;
    if numel(unique(xValues)) > 1
        lineFit = polyfit(xValues,yValues,1);
        plotX = linspace(min(xValues),max(xValues),100);
        plot(axLink,plotX,polyval(lineFit,plotX), ...
            '-','Color',[0.40 0.40 0.40],'LineWidth',1.5);
    end
    xlabel(axLink,'PPS Post-training - Baseline (m)');
    ylabel(axLink,'SpeedPad benefit (hits)');
    title(axLink,sprintf('Spearman rho = %.2f, exact p = %.3f', ...
        spearmanRho,spearmanP));
    local_styleAxes(axLink);

    exportgraphics(figIntegrated,fullfile(outputFolder, ...
        'WP4_Integrated_PPS_SpeedPad_Figure.png'),'Resolution',300);
    savefig(figIntegrated,fullfile(outputFolder, ...
        'WP4_Integrated_PPS_SpeedPad_Figure.fig'));
end

%% ======================== COMMAND-WINDOW REPORT ======================

fprintf('\n============================================================\n');
fprintf('WP4 INTEGRATED PRELIMINARY RESULTS\n');
fprintf('============================================================\n');
fprintf('Participants included: %d\n',nParticipants);
fprintf('Included IDs:\n');
disp(IntegratedData.ParticipantCode);

fprintf('\nParticipant-level PPS:\n');
fprintf('  Baseline:      M = %.4f, SD = %.4f m\n', ...
    mean(IntegratedData.BaselinePPSBorder_m), ...
    std(IntegratedData.BaselinePPSBorder_m));
fprintf('  Post-training: M = %.4f, SD = %.4f m\n', ...
    mean(IntegratedData.PostTrainingPPSBorder_m), ...
    std(IntegratedData.PostTrainingPPSBorder_m));
fprintf('  Post-training - baseline: M = %.4f, SD = %.4f m\n', ...
    ppsMeanChange,std(IntegratedData.PPSChange_m));
fprintf('  Exact two-sided sign-flip p = %.5f\n',ppsPtwo);
fprintf('  Directional p (Post-training > baseline) = %.5f\n', ...
    ppsPgreater);

fprintf('\nSpeedPad Time-by-Condition contrast:\n');
fprintf('  Pull Training change: M = %.2f, SD = %.2f hits\n', ...
    mean(IntegratedData.PullTrainingChange), ...
    std(IntegratedData.PullTrainingChange));
fprintf('  Control change: M = %.2f, SD = %.2f hits\n', ...
    mean(IntegratedData.ControlChange),std(IntegratedData.ControlChange));
fprintf('  Pull Training-specific benefit: M = %.2f, SD = %.2f hits\n', ...
    speedMeanBenefit,std(IntegratedData.SpeedPadBenefit_hits));
fprintf('  Exact two-sided sign-flip p = %.5f\n',speedPtwo);

fprintf('\nPPS change versus SpeedPad benefit:\n');
fprintf('  Exact Spearman rho = %.4f, p = %.5f\n', ...
    spearmanRho,spearmanP);

if any(IntegratedData.BaselineFitR2Below_0_80 | ...
       IntegratedData.PostTrainingFitR2Below_0_80)
    fprintf(['\nFIT-QUALITY NOTE: One or more participant-level sigmoid ' ...
        'fits had R2 < %.2f.\nThese participants were retained because ' ...
        'the threshold is reported as a transparency flag, not as a ' ...
        'post-training exclusion rule.\n'],fitQualityFlagThreshold);
    disp(IntegratedData( ...
        IntegratedData.BaselineFitR2Below_0_80 | ...
        IntegratedData.PostTrainingFitR2Below_0_80, ...
        {'ParticipantCode','BaselineSigmoidR2', ...
         'PostTrainingSigmoidR2'}));
end

fprintf(['\nAll findings are preliminary. The PPS-SpeedPad association ' ...
    'must be interpreted as exploratory at the current sample size.\n']);
fprintf('Outputs saved to:\n%s\n',outputFolder);
fprintf('============================================================\n');

% Save only numerical and tabular outputs. Excluding graphics handles avoids
% the large-file warning that occurs when the complete workspace is saved.
save(fullfile(outputFolder,'WP4_Integrated_Analysis_Workspace.mat'), ...
    'IntegratedData','PPSSummaryAll','SpeedPadAll','MainTests', ...
    'AssociationResults','distances','ppsMeanChange','ppsPtwo', ...
    'ppsPgreater','speedMeanBenefit','speedPtwo','spearmanRho', ...
    'spearmanP','fitQualityFlagThreshold');

%% ======================== LOCAL FUNCTIONS ============================

function participantID = local_numericParticipantID(participantCode)
    digitsOnly = regexprep(participantCode,'\D','');
    participantID = str2double(digitsOnly);
    if ~isfinite(participantID)
        error('Could not extract a numeric ID from %s.',participantCode);
    end
end

function filePath = local_resolvePPSFile(folderPath,participantCode,phase)
    pattern = [participantCode phase '*.csv'];
    candidates = dir(fullfile(folderPath,pattern));
    validFiles = {};

    for i = 1:numel(candidates)
        candidatePath = fullfile(candidates(i).folder,candidates(i).name);
        if local_isBehaviouralPPSFile(candidatePath)
            validFiles{end+1,1} = candidatePath; %#ok<AGROW>
        end
    end

    if isempty(validFiles)
        error(['No behavioural PPS file was found for %s%s in:\n%s\n' ...
            'The required columns are Vibrate_distance, RT, and Monster.'], ...
            participantCode,phase,folderPath);
    end
    if numel(validFiles) > 1
        error(['More than one behavioural PPS file was found for %s%s:\n%s\n' ...
            'Keep only one valid behavioural file per phase.'], ...
            participantCode,phase,strjoin(validFiles,newline));
    end
    filePath = validFiles{1};
end

function isValid = local_isBehaviouralPPSFile(filePath)
    fileID = fopen(filePath,'r');
    if fileID < 0
        isValid = false;
        return
    end
    cleanUp = onCleanup(@() fclose(fileID)); %#ok<NASGU>
    header = lower(fgetl(fileID));
    isValid = contains(header,'vibrate_distance') && ...
              contains(header,'rt') && contains(header,'monster');
end

function [meanRT,nValid] = local_conditionMeans( ...
    filePath,distances,minimumRT,maximumRT,outlierCutoff)

    T = readtable(filePath,'VariableNamingRule','preserve');
    trialDistance = local_getNumericColumn(T,'Vibrate_distance');
    rt = local_getNumericColumn(T,'RT');
    monster = local_getNumericColumn(T,'Monster');

    meanRT = nan(1,numel(distances));
    nValid = zeros(1,numel(distances));

    for d = 1:numel(distances)
        use = monster == 0 & ...
              abs(trialDistance-distances(d)) < 1e-9 & ...
              isfinite(rt) & rt >= minimumRT & rt <= maximumRT;
        trials = rt(use);

        if isempty(trials)
            continue
        end

        trialMean = mean(trials);
        trialSD = std(trials);
        if isfinite(trialSD) && trialSD > 0
            zScore = abs((trials-trialMean)/trialSD);
            trials = trials(zScore < outlierCutoff);
        end

        nValid(d) = sum(isfinite(trials));
        if nValid(d) > 0
            meanRT(d) = mean(trials,'omitnan');
        end
    end

    if any(~isfinite(meanRT))
        error('Incomplete PPS distance means in %s.',filePath);
    end
end

function [border,rSquared] = local_fitPPSBorder( ...
    distances,meanRT,eqsig,fitOptions)

    border = NaN;
    rSquared = NaN;

    if any(~isfinite(meanRT)) || max(meanRT) <= min(meanRT)
        return
    end

    normalisedRT = (meanRT-min(meanRT))/(max(meanRT)-min(meanRT));

    try
        [fitObject,goodness] = fit( ...
            distances,normalisedRT,eqsig,fitOptions);
        border = fitObject.x0;
        rSquared = goodness.rsquare;
    catch
        border = NaN;
        rSquared = NaN;
    end
end

function values = local_getNumericColumn(T,requestedName)
    names = T.Properties.VariableNames;
    normalisedNames = cell(size(names));
    for i = 1:numel(names)
        normalisedNames{i} = lower(regexprep(strtrim(names{i}),'[^a-zA-Z0-9]',''));
    end
    requested = lower(regexprep(strtrim(requestedName),'[^a-zA-Z0-9]',''));
    index = find(strcmp(normalisedNames,requested),1);
    if isempty(index)
        error('Required column "%s" was not found.',requestedName);
    end
    raw = T.(names{index});
    if isnumeric(raw) || islogical(raw)
        values = double(raw(:));
    else
        values = str2double(string(raw(:)));
    end
end

function numericValues = local_extractCounts( ...
    rawColumn,requireComplete,variableLabel)

    numericValues = nan(numel(rawColumn),1);
    for rowNumber = 1:numel(rawColumn)
        currentValue = rawColumn{rowNumber};
        if isnumeric(currentValue) && isscalar(currentValue)
            if isfinite(currentValue)
                numericValues(rowNumber) = double(currentValue);
            end
            continue
        end

        currentText = strtrim(string(currentValue));
        if ismissing(currentText) || strlength(currentText) == 0
            continue
        end

        firstNumber = regexp(char(currentText), ...
            '[-+]?(?:\d+\.?\d*|\.\d+)','match','once');
        if ~isempty(firstNumber)
            numericValues(rowNumber) = str2double(firstNumber);
        end
    end

    if requireComplete && any(~isfinite(numericValues))
        badRows = find(~isfinite(numericValues));
        error('Could not extract %s in row(s): %s.', ...
            variableLabel,mat2str(badRows'));
    end
end

function [observedMean,pTwoSided,pGreater,pLess, ...
    nPermutations,method] = local_signFlip(differences)

    differences = differences(:);
    differences = differences(isfinite(differences));
    sampleSize = numel(differences);
    if sampleSize < 1
        error('The sign-flip test requires at least one contrast.');
    end

    observedMean = mean(differences);
    tolerance = 1e-12*max(1,abs(observedMean));

    if sampleSize <= 20
        nPermutations = 2^sampleSize;
        method = "Exact enumeration";
        permutationStatistics = nan(nPermutations,1);
        for permutation = 0:(nPermutations-1)
            bits = bitget(uint64(permutation),1:sampleSize);
            signs = 2*double(bits(:))-1;
            permutationStatistics(permutation+1) = ...
                mean(differences.*signs);
        end
        pTwoSided = mean(abs(permutationStatistics) >= ...
            abs(observedMean)-tolerance);
        pGreater = mean(permutationStatistics >= observedMean-tolerance);
        pLess = mean(permutationStatistics <= observedMean+tolerance);
    else
        nPermutations = 100000;
        method = "Monte Carlo approximation";
        randomSigns = 2*(rand(nPermutations,sampleSize) >= 0.5)-1;
        permutationStatistics = mean(randomSigns.*differences',2);
        pTwoSided = (sum(abs(permutationStatistics) >= ...
            abs(observedMean)-tolerance)+1)/(nPermutations+1);
        pGreater = (sum(permutationStatistics >= ...
            observedMean-tolerance)+1)/(nPermutations+1);
        pLess = (sum(permutationStatistics <= ...
            observedMean+tolerance)+1)/(nPermutations+1);
    end
end

function [rho,pValue,nPermutations,method] = ...
    local_exactSpearman(x,y)

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    sampleSize = numel(x);
    if sampleSize < 3
        error('Spearman analysis requires at least three complete pairs.');
    end

    rankX = local_tiedRanks(x);
    rankY = local_tiedRanks(y);
    rho = local_pearson(rankX,rankY);
    tolerance = 1e-12*max(1,abs(rho));

    if sampleSize <= 9
        permutationIndex = perms(1:sampleSize);
        nPermutations = size(permutationIndex,1);
        method = "Exact enumeration";
        permutationRho = nan(nPermutations,1);
        for permutation = 1:nPermutations
            permutationRho(permutation) = local_pearson( ...
                rankX,rankY(permutationIndex(permutation,:)));
        end
        pValue = mean(abs(permutationRho) >= abs(rho)-tolerance);
    else
        nPermutations = 100000;
        method = "Monte Carlo approximation";
        extremeCount = 0;
        for permutation = 1:nPermutations
            permutedY = rankY(randperm(sampleSize));
            permutedRho = local_pearson(rankX,permutedY);
            extremeCount = extremeCount + ...
                (abs(permutedRho) >= abs(rho)-tolerance);
        end
        pValue = (extremeCount+1)/(nPermutations+1);
    end
end

function ranks = local_tiedRanks(values)
    values = values(:);
    [sortedValues,order] = sort(values);
    sortedRanks = nan(size(values));
    first = 1;
    while first <= numel(values)
        last = first;
        while last < numel(values) && ...
                sortedValues(last+1) == sortedValues(first)
            last = last+1;
        end
        sortedRanks(first:last) = mean(first:last);
        first = last+1;
    end
    ranks = nan(size(values));
    ranks(order) = sortedRanks;
end

function value = local_pearson(x,y)
    x = x(:)-mean(x);
    y = y(:)-mean(y);
    denominator = sqrt(sum(x.^2)*sum(y.^2));
    if denominator <= 0
        value = NaN;
    else
        value = sum(x.*y)/denominator;
    end
end

function local_groupMeanSEM(targetAxes,xPosition,values,color)
    values = values(isfinite(values));
    groupMean = mean(values);
    groupSEM = std(values)/sqrt(numel(values));
    errorbar(targetAxes,xPosition,groupMean,groupSEM,'o', ...
        'Color',color,'MarkerFaceColor',color,'MarkerSize',8, ...
        'LineWidth',2,'CapSize',8,'HandleVisibility','off');
end

function local_styleAxes(targetAxes)
    box(targetAxes,'off');
    pbaspect(targetAxes,[1 1 1]);
    set(targetAxes,'Color','w','FontName','Arial','FontSize',12, ...
        'LineWidth',2,'TickDir','out','TickLength',[0.025 0.025], ...
        'XColor','k','YColor','k','Layer','top');
end

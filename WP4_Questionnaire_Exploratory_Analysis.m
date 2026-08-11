%% WP4 QUESTIONNAIRE EXPLORATORY ANALYSIS
% Standalone participant-level analysis linking questionnaire total scores
% with the two principal WP4 change outcomes:
%   1. PPS change (Post-training minus baseline, in metres)
%   2. SpeedPad benefit beyond control
%      (Pull Training Second Assessment - First Assessment) minus
%      (Control Second Assessment - First Assessment)
%
% IMPORTANT:
%   - The script performs the PPS cleaning and predefined baseline QC
%     internally, so no other WP4 analysis script must be run first.
%   - Associations are exploratory because the current sample is small.
%   - Exact Spearman permutation tests are used when N <= 9. For larger
%     future samples, a reproducible Monte Carlo permutation test is used.
%   - Benjamini-Hochberg FDR correction is applied across all questionnaire
%     association tests.
%   - Statistics and Machine Learning Toolbox is not required.
%   - Subject-level PPS sigmoid fitting uses the same Curve Fitting Toolbox
%     functions already used in the WP3/WP4 PPS analyses.

clear; close all; clc;

dbstop if error
% dbstop if warning

%% ======================== USER SETTINGS ==============================

dataFolder = '';

questionnaireFile = fullfile(dataFolder,'WP4_Data Base.xlsx');
questionnaireSheet = 'ALL';

speedPadDataFile = fullfile(dataFolder,'speedpabdata.xlsx');

% Candidate participant prefixes and PPS source-file patterns.
participantPrefixes = { ...
    ''};
ppsFilePatterns = {'pre*.csv','post*.csv'};
% User-facing labels. File patterns remain pre/post because they correspond
% to the existing CSV filenames.
ppsConditionNames = {'Baseline','Post-training'};
distances = [0.25; 0.75; 1.25; 1.75; 2.25];

% Same cleaning and baseline-QC settings used in the WP4 PPS script.
centralonly = 0;      % 0 = central stimuli only
outlierCutoff = 3;    % within-person/condition/distance z threshold
minimumRT = 0.1;      % seconds
maximumRT = 1.0;      % seconds
minimumBaselineR2 = 0.10;

makeFigures = true;
alpha = 0.05;

%% ================= PPS CLEANING AND BASELINE QC ======================

name = participantPrefixes(:);
nCandidateParticipants = numel(name);
nConditions = numel(ppsFilePatterns);
nDistances = numel(distances);

if nConditions ~= 2
    error(['Exactly two PPS assessments (baseline and Post-training) ' ...
        'are required.']);
end

rtAll = nan(nCandidateParticipants,nConditions,nDistances);
baselineSlope = nan(nCandidateParticipants,1);
baselineR2 = nan(nCandidateParticipants,1);

for participantNumber = 1:nCandidateParticipants
    for conditionNumber = 1:nConditions
        files = dir(fullfile(dataFolder,[ ...
            name{participantNumber} ...
            ppsFilePatterns{conditionNumber}]));

        if isempty(files)
            error('Missing PPS file for %s, condition %s.', ...
                name{participantNumber}, ...
                ppsConditionNames{conditionNumber});
        end

        if numel(files) > 1
            error(['More than one PPS file was found for %s, ' ...
                'condition %s.'],name{participantNumber}, ...
                ppsConditionNames{conditionNumber});
        end

        sourceTable = readtable(fullfile( ...
            files(1).folder,files(1).name), ...
            'VariableNamingRule','preserve');

        requiredVariables = {'Vibrate_distance','RT','Monster'};

        if ~all(ismember(requiredVariables, ...
                sourceTable.Properties.VariableNames))
            error(['Required PPS variables are missing from %s. ' ...
                'Expected Vibrate_distance, RT, and Monster.'], ...
                files(1).name);
        end

        currentData = [ ...
            sourceTable.Vibrate_distance, ...
            sourceTable.RT, ...
            sourceTable.Monster];

        if centralonly == 0
            currentData = currentData(currentData(:,3) == 0,:);
        elseif centralonly == 2
            currentData = currentData(currentData(:,3) ~= 0,:);
        end

        currentData(currentData(:,1) == -1,2) = NaN;
        currentData(currentData(:,1) == 0,2) = NaN;
        currentData(currentData(:,2) < minimumRT,2) = NaN;
        currentData(currentData(:,2) > maximumRT,2) = NaN;
        currentData = currentData(isfinite(currentData(:,2)),:);

        for distanceNumber = 1:nDistances
            trials = currentData( ...
                currentData(:,1) == distances(distanceNumber),2);

            if isempty(trials)
                continue
            end

            trialMean = mean(trials,'omitnan');
            trialSD = std(trials,'omitnan');

            if isfinite(trialSD) && trialSD > 0
                trialZ = abs((trials-trialMean)/trialSD);
                trials = trials(trialZ < outlierCutoff);
            end

            rtAll(participantNumber,conditionNumber,distanceNumber) = ...
                mean(trials,'omitnan');
        end
    end

    baselineRT = squeeze(rtAll(participantNumber,1,:));

    if all(isfinite(baselineRT)) && max(baselineRT) > min(baselineRT)
        linearCoefficients = polyfit(distances,baselineRT,1);
        predictedRT = polyval(linearCoefficients,distances);
        residualSS = sum((baselineRT-predictedRT).^2);
        totalSS = sum((baselineRT-mean(baselineRT)).^2);

        if totalSS > 0
            baselineSlope(participantNumber) = linearCoefficients(1);
            baselineR2(participantNumber) = 1-residualSS/totalSS;
        end
    end
end

lowOrInvalidR2 = baselineR2 < minimumBaselineR2 | ...
    ~isfinite(baselineR2);
nonPositiveSlope = baselineSlope <= 0 | ~isfinite(baselineSlope);
excludedByPPSQC = lowOrInvalidR2 | nonPositiveSlope;

qcTable = table( ...
    string(name),baselineSlope,baselineR2, ...
    lowOrInvalidR2,nonPositiveSlope,excludedByPPSQC, ...
    'VariableNames',{ ...
    'Subject','BaselineSlope','BaselineR2', ...
    'LowOrInvalid_BaselineR2','NonPositive_BaselineSlope', ...
    'Excluded_BaselineQC'});

writetable(qcTable,fullfile( ...
    dataFolder,'WP4_Questionnaire_PPS_QC_summary.csv'));

fprintf('\nPPS baseline QC completed. Excluded participant(s):\n');
disp(qcTable(qcTable.Excluded_BaselineQC,:));

name(excludedByPPSQC) = [];
rtAll(excludedByPPSQC,:,:) = [];

rtPre = squeeze(rtAll(:,1,:));
rtPost = squeeze(rtAll(:,2,:));
nPPSParticipants = numel(name);

if nPPSParticipants < 3
    error('Fewer than three participants remained after PPS QC.');
end

participantID = nan(nPPSParticipants,1);

for participantNumber = 1:nPPSParticipants
    idText = regexp(name{participantNumber},'^\d+','match','once');

    if isempty(idText)
        error('Could not extract a numeric ID from %s.', ...
            name{participantNumber});
    end

    participantID(participantNumber) = str2double(idText);
end

if numel(unique(participantID)) ~= nPPSParticipants
    error('Duplicate participant IDs were detected in the PPS workspace.');
end

%% ================= SUBJECT-LEVEL PPS BORDERS =========================

eqsig = fittype( ...
    '1/(1+exp(-(x-x0)/b))', ...
    'independent','x', ...
    'coefficients',{'b','x0'});

fitOptions = fitoptions(eqsig);
fitOptions.Lower = [eps 0.25];
fitOptions.Upper = [Inf 2.25];
fitOptions.StartPoint = [0.5 1.25];

ppsPreBorder = nan(nPPSParticipants,1);
ppsPostBorder = nan(nPPSParticipants,1);
ppsPreFitR2 = nan(nPPSParticipants,1);
ppsPostFitR2 = nan(nPPSParticipants,1);

for participantNumber = 1:nPPSParticipants
    [ppsPreBorder(participantNumber), ...
     ppsPreFitR2(participantNumber)] = fit_subject_pps_border( ...
        distances,rtPre(participantNumber,:)',eqsig,fitOptions);

    [ppsPostBorder(participantNumber), ...
     ppsPostFitR2(participantNumber)] = fit_subject_pps_border( ...
        distances,rtPost(participantNumber,:)',eqsig,fitOptions);
end

ppsChangeMetres = ppsPostBorder-ppsPreBorder;
ppsPercentChange = 100*ppsChangeMetres./ppsPreBorder;

ppsParticipantTable = table( ...
    participantID,string(name), ...
    ppsPreBorder,ppsPostBorder,ppsChangeMetres,ppsPercentChange, ...
    ppsPreFitR2,ppsPostFitR2, ...
    'VariableNames',{ ...
    'ParticipantID','PPS_FilePrefix', ...
    'PPS_Baseline_m','PPS_PostTraining_m','PPS_Change_m', ...
    'PPS_PercentChange','PPS_Baseline_FitR2', ...
    'PPS_PostTraining_FitR2'});

%% ======================== LOAD SPEEDPAD ==============================

if ~isfile(speedPadDataFile)
    error('The raw SpeedPad workbook was not found:\n%s', ...
        speedPadDataFile);
end

rawSpeedPad = readcell(speedPadDataFile);

if size(rawSpeedPad,1) < 2 || size(rawSpeedPad,2) < 14
    error(['The SpeedPad workbook does not match the expected layout ' ...
        'in columns A:N.']);
end

speedPadIDCandidates = numeric_cell_column( ...
    rawSpeedPad(2:end,1),'SpeedPad participant ID');
speedPadRowsMask = isfinite(speedPadIDCandidates);
speedPadRows = rawSpeedPad(1+find(speedPadRowsMask),:);
speedPadID = speedPadIDCandidates(speedPadRowsMask);

pullTrainingFirstAssessment = speedpad_count_column( ...
    speedPadRows(:,3),'Pull Training First Assessment Hits');
pullTrainingSecondAssessment = speedpad_count_column( ...
    speedPadRows(:,6),'Pull Training Second Assessment Hits');
controlFirstAssessment = speedpad_count_column( ...
    speedPadRows(:,10),'Control First Assessment Hits');
controlSecondAssessment = speedpad_count_column( ...
    speedPadRows(:,13),'Control Second Assessment Hits');

speedPadBenefit = ...
    (pullTrainingSecondAssessment-pullTrainingFirstAssessment) - ...
    (controlSecondAssessment-controlFirstAssessment);

speedPadTable = table( ...
    speedPadID,speedPadBenefit, ...
    'VariableNames',{'ParticipantID','SpeedPad_Benefit_Hits'});

if any(~isfinite(speedPadTable.SpeedPad_Benefit_Hits))
    error('Missing or invalid SpeedPad hit values were detected.');
end

if numel(unique(speedPadTable.ParticipantID)) ~= height(speedPadTable)
    error('Duplicate participant IDs were found in the SpeedPad file.');
end

%% ==================== LOAD QUESTIONNAIRE TOTALS ======================

if ~isfile(questionnaireFile)
    error('The questionnaire workbook was not found:\n%s', ...
        questionnaireFile);
end


rawQuestionnaires = readcell( ...
    questionnaireFile,'Sheet',questionnaireSheet);

% Expected workbook structure (two header rows):
%   A  (1)  = Subject Number
%   AC (29) = Total ASRS
%   AZ (52) = ATTC Total
%   BQ (69) = TEQ Total
%   CH (86) = Final Total SSQ score
requiredQuestionnaireColumns = [1 29 52 69 86];

if size(rawQuestionnaires,1) < 3 || ...
        size(rawQuestionnaires,2) < max(requiredQuestionnaireColumns)
    error(['The questionnaire workbook does not match the expected ' ...
        'ALL-sheet layout with two header rows and columns A:CH.']);
end

questionnaireRows = rawQuestionnaires(3:end,:);

questionnaireID = numeric_cell_column( ...
    questionnaireRows(:,1),'Subject Number');
asrsTotal = numeric_cell_column( ...
    questionnaireRows(:,29),'Total ASRS');
attcTotal = numeric_cell_column( ...
    questionnaireRows(:,52),'ATTC Total');
teqTotal = numeric_cell_column( ...
    questionnaireRows(:,69),'TEQ Total');
ssqTotal = numeric_cell_column( ...
    questionnaireRows(:,86),'Final Total SSQ score');

populatedRows = isfinite(questionnaireID);

questionnaireTable = table( ...
    questionnaireID(populatedRows), ...
    asrsTotal(populatedRows), ...
    attcTotal(populatedRows), ...
    teqTotal(populatedRows), ...
    ssqTotal(populatedRows), ...
    'VariableNames',{ ...
    'ParticipantID','ASRS_Total','ATTC_Total','TEQ_Total','SSQ_Total'});

if numel(unique(questionnaireTable.ParticipantID)) ~= ...
        height(questionnaireTable)
    error('Duplicate participant IDs were found in the questionnaire file.');
end

%% ======================== MERGE BY ID ================================

analysisData = innerjoin( ...
    ppsParticipantTable,speedPadTable,'Keys','ParticipantID');

analysisData = innerjoin( ...
    analysisData,questionnaireTable,'Keys','ParticipantID');

analysisData = sortrows(analysisData,'ParticipantID');

if height(analysisData) ~= nPPSParticipants
    missingIDs = setdiff(participantID,analysisData.ParticipantID);
    error(['Not all PPS-QC-approved participants could be matched across ' ...
        'PPS, SpeedPad, and questionnaire data. Missing IDs: %s'], ...
        mat2str(missingIDs'));
end

fprintf('\nIncluded participant IDs after PPS QC and data matching:\n');
disp(analysisData.ParticipantID');

%% ================= DESCRIPTIVE STATISTICS ============================

descriptiveVariables = { ...
    'ASRS_Total','ATTC_Total','TEQ_Total','SSQ_Total', ...
    'PPS_Baseline_m','PPS_PostTraining_m','PPS_Change_m', ...
    'SpeedPad_Benefit_Hits'};

descriptiveLabels = { ...
    'ASRS total','ATTC total','TEQ total','SSQ total', ...
    'PPS baseline (m)','PPS Post-training (m)','PPS change (m)', ...
    'SpeedPad benefit (hits)'};

nDescriptive = numel(descriptiveVariables);
descriptiveN = nan(nDescriptive,1);
descriptiveMean = nan(nDescriptive,1);
descriptiveSD = nan(nDescriptive,1);
descriptiveMedian = nan(nDescriptive,1);
descriptiveMinimum = nan(nDescriptive,1);
descriptiveMaximum = nan(nDescriptive,1);

for variableNumber = 1:nDescriptive
    currentValues = analysisData.( ...
        descriptiveVariables{variableNumber});
    currentValues = currentValues(isfinite(currentValues));

    descriptiveN(variableNumber) = numel(currentValues);
    descriptiveMean(variableNumber) = mean(currentValues);
    descriptiveSD(variableNumber) = std(currentValues);
    descriptiveMedian(variableNumber) = median(currentValues);
    descriptiveMinimum(variableNumber) = min(currentValues);
    descriptiveMaximum(variableNumber) = max(currentValues);
end

descriptiveTable = table( ...
    string(descriptiveLabels(:)), ...
    descriptiveN,descriptiveMean,descriptiveSD,descriptiveMedian, ...
    descriptiveMinimum,descriptiveMaximum, ...
    'VariableNames',{ ...
    'Variable','N','Mean','SD','Median','Minimum','Maximum'});

%% ============= EXACT SPEARMAN PERMUTATION TESTS =====================

questionnaireVariables = { ...
    'ASRS_Total','ATTC_Total','TEQ_Total','SSQ_Total'};

questionnaireLabels = { ...
    'ASRS','ATTC','TEQ','SSQ'};

outcomeVariables = { ...
    'PPS_Change_m','SpeedPad_Benefit_Hits'};

outcomeLabels = { ...
    'PPS change: Post-training - baseline (m)', ...
    'SpeedPad benefit (hits)'};

% Shorter labels are used on the figure axes for readability.
outcomeAxisLabels = { ...
    'PPS change (m)', ...
    'SpeedPad benefit (hits)'};

nTests = numel(questionnaireVariables)*numel(outcomeVariables);
resultQuestionnaire = strings(nTests,1);
resultOutcome = strings(nTests,1);
resultN = nan(nTests,1);
resultRho = nan(nTests,1);
resultP = nan(nTests,1);
resultPermutations = nan(nTests,1);
resultMethod = strings(nTests,1);

testNumber = 0;

for questionnaireNumber = 1:numel(questionnaireVariables)
    for outcomeNumber = 1:numel(outcomeVariables)
        testNumber = testNumber+1;

        x = analysisData.( ...
            questionnaireVariables{questionnaireNumber});
        y = analysisData.( ...
            outcomeVariables{outcomeNumber});

        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);

        [rho,pValue,nPermutations,method] = ...
            spearman_permutation_base(x,y);

        resultQuestionnaire(testNumber) = ...
            questionnaireLabels{questionnaireNumber};
        resultOutcome(testNumber) = outcomeLabels{outcomeNumber};
        resultN(testNumber) = numel(x);
        resultRho(testNumber) = rho;
        resultP(testNumber) = pValue;
        resultPermutations(testNumber) = nPermutations;
        resultMethod(testNumber) = method;
    end
end

fdrP = benjamini_hochberg_base(resultP);
rejectFDR = fdrP < alpha;

correlationTable = table( ...
    resultQuestionnaire,resultOutcome,resultN,resultRho,resultP,fdrP, ...
    rejectFDR,resultPermutations,resultMethod, ...
    'VariableNames',{ ...
    'Questionnaire','Outcome','N','SpearmanRho','ExactPermutationP', ...
    'FDR_P','RejectFDR_0_05','NumberOfPermutations','PermutationMethod'});

%% ======================== DISPLAY RESULTS ============================

fprintf('\n============================================================\n');
fprintf('WP4 QUESTIONNAIRE EXPLORATORY RESULTS\n');
fprintf('============================================================\n');
fprintf('Participants included: %d\n',height(analysisData));
fprintf(['Participant 1 and any future PPS-QC exclusions are removed ' ...
    'automatically through the internal baseline-QC procedure.\n\n']);

disp(descriptiveTable);

fprintf('\nExploratory exact Spearman associations:\n');
disp(correlationTable);

fprintf(['\nInterpret all associations as exploratory. FDR correction ' ...
    'was applied across %d tests.\n'],nTests);
fprintf('============================================================\n');

%% ======================== SAVE OUTPUTS ===============================

writetable(analysisData,fullfile( ...
    dataFolder,'WP4_Questionnaire_participant_data.csv'));

writetable(descriptiveTable,fullfile( ...
    dataFolder,'WP4_Questionnaire_descriptive_summary.csv'));

writetable(correlationTable,fullfile( ...
    dataFolder,'WP4_Questionnaire_exact_Spearman_results.csv'));

save(fullfile(dataFolder,'WP4_Questionnaire_results.mat'), ...
    'analysisData','descriptiveTable','correlationTable', ...
    'questionnaireVariables','questionnaireLabels','outcomeVariables', ...
    'outcomeLabels','outcomeAxisLabels','alpha');

%% ======================== FIGURES ====================================

if makeFigures
    figureQuestionnaires = figure(1); clf(figureQuestionnaires);
    set(figureQuestionnaires,'Color','w','Position',[80 80 900 1040]);

    plotLayout = tiledlayout( ...
        figureQuestionnaires,4,2, ...
        'TileSpacing','compact','Padding','compact');

    pointColour = [0.18 0.28 0.70];
    lineColour = [0.45 0.45 0.45];
    testNumber = 0;

    for questionnaireNumber = 1:numel(questionnaireVariables)
        for outcomeNumber = 1:numel(outcomeVariables)
            testNumber = testNumber+1;
            targetAxes = nexttile(plotLayout);
            hold(targetAxes,'on');

            x = analysisData.( ...
                questionnaireVariables{questionnaireNumber});
            y = analysisData.( ...
                outcomeVariables{outcomeNumber});
            valid = isfinite(x) & isfinite(y);
            x = x(valid);
            y = y(valid);

            scatter(targetAxes,x,y,52, ...
                'MarkerFaceColor','w', ...
                'MarkerEdgeColor',pointColour, ...
                'LineWidth',1.4);

            % Descriptive linear trend only; inference is based on the
            % exact Spearman permutation test reported in the title.
            if numel(unique(x)) > 1
                trendCoefficients = polyfit(x,y,1);
                trendX = linspace(min(x),max(x),100);
                trendY = polyval(trendCoefficients,trendX);
                plot(targetAxes,trendX,trendY,'-', ...
                    'Color',lineColour,'LineWidth',1.3, ...
                    'HandleVisibility','off');
            end

            xlabel(targetAxes,questionnaireLabels{questionnaireNumber});
            ylabel(targetAxes,outcomeAxisLabels{outcomeNumber});
            title(targetAxes,sprintf( ...
                '\\rho = %.2f, exact p = %.3f, FDR p = %.3f', ...
                resultRho(testNumber),resultP(testNumber), ...
                fdrP(testNumber)), ...
                'FontWeight','normal','FontSize',10);

            style_questionnaire_axes(targetAxes);
        end
    end

    title(plotLayout, ...
        'WP4 exploratory questionnaire associations', ...
        'FontName','Arial','FontSize',14,'FontWeight','normal');

    exportgraphics(figureQuestionnaires,fullfile( ...
        dataFolder,'WP4_Questionnaire_associations.png'), ...
        'Resolution',300);
end

fprintf('\nAll questionnaire outputs were saved to:\n%s\n',dataFolder);

%% ======================== LOCAL FUNCTIONS ============================

function [ppsBorder,fitR2] = fit_subject_pps_border( ...
    distances,meanRT,eqsig,fitOptions)

    distances = distances(:);
    meanRT = meanRT(:);

    if any(~isfinite(meanRT)) || max(meanRT) <= min(meanRT)
        ppsBorder = NaN;
        fitR2 = NaN;
        return
    end

    normalisedRT = ...
        (meanRT-min(meanRT))/(max(meanRT)-min(meanRT));

    try
        [fitObject,goodness] = fit( ...
            distances,normalisedRT,eqsig,fitOptions);
        ppsBorder = fitObject.x0;
        fitR2 = goodness.rsquare;
    catch
        ppsBorder = NaN;
        fitR2 = NaN;
    end
end

function numericValues = numeric_cell_column(rawColumn,columnLabel)

    numericValues = nan(numel(rawColumn),1);

    for rowNumber = 1:numel(rawColumn)
        currentValue = rawColumn{rowNumber};

        if isnumeric(currentValue) && isscalar(currentValue)
            numericValues(rowNumber) = double(currentValue);
        elseif islogical(currentValue) && isscalar(currentValue)
            numericValues(rowNumber) = double(currentValue);
        elseif ischar(currentValue) || ...
                (isstring(currentValue) && isscalar(currentValue))
            parsedValue = str2double(strtrim(string(currentValue)));

            if isfinite(parsedValue)
                numericValues(rowNumber) = parsedValue;
            end
        elseif ~isempty(currentValue)
            warning('Non-numeric value ignored in %s, data row %d.', ...
                columnLabel,rowNumber);
        end
    end
end

function numericValues = speedpad_count_column(rawColumn,columnLabel)

    % Accept both numeric cells and strings such as "122 (0.0%)".
    numericValues = nan(numel(rawColumn),1);

    for rowNumber = 1:numel(rawColumn)
        currentValue = rawColumn{rowNumber};

        if isnumeric(currentValue) && isscalar(currentValue)
            numericValues(rowNumber) = double(currentValue);
        elseif islogical(currentValue) && isscalar(currentValue)
            numericValues(rowNumber) = double(currentValue);
        elseif ischar(currentValue) || ...
                (isstring(currentValue) && isscalar(currentValue))
            numberToken = regexp(char(string(currentValue)), ...
                '[-+]?\d*\.?\d+','match','once');

            if ~isempty(numberToken)
                numericValues(rowNumber) = str2double(numberToken);
            end
        end

        if ~isfinite(numericValues(rowNumber))
            error('Invalid value in %s, data row %d.', ...
                columnLabel,rowNumber);
        end
    end
end

function [rho,pValue,nPermutations,method] = ...
    spearman_permutation_base(x,y)

    x = x(:);
    y = y(:);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    sampleSize = numel(x);

    if sampleSize < 3
        error('At least three complete observations are required.');
    end

    rankX = average_ranks_base(x);
    rankY = average_ranks_base(y);
    rho = pearson_base(rankX,rankY);

    if ~isfinite(rho)
        pValue = NaN;
        nPermutations = 0;
        method = "Not calculated: constant variable";
        return
    end

    tolerance = 1e-12*max(1,abs(rho));

    if sampleSize <= 9
        permutationIndices = perms(1:sampleSize);
        nPermutations = size(permutationIndices,1);
        extremeCount = 0;

        for permutationNumber = 1:nPermutations
            permutedRho = pearson_base( ...
                rankX,rankY(permutationIndices(permutationNumber,:)));
            extremeCount = extremeCount + ...
                (abs(permutedRho) >= abs(rho)-tolerance);
        end

        pValue = extremeCount/nPermutations;
        method = "Exact enumeration";
    else
        nPermutations = 100000;
        extremeCount = 0;
        rng(20260804,'twister');

        for permutationNumber = 1:nPermutations
            permutedRho = pearson_base( ...
                rankX,rankY(randperm(sampleSize)));
            extremeCount = extremeCount + ...
                (abs(permutedRho) >= abs(rho)-tolerance);
        end

        pValue = (extremeCount+1)/(nPermutations+1);
        method = "Monte Carlo permutation";
    end

    pValue = min(max(pValue,0),1);
end

function ranks = average_ranks_base(values)

    values = values(:);
    [sortedValues,sortOrder] = sort(values);
    ranks = nan(size(values));
    firstIndex = 1;

    while firstIndex <= numel(values)
        lastIndex = firstIndex;

        while lastIndex < numel(values) && ...
                sortedValues(lastIndex+1) == sortedValues(firstIndex)
            lastIndex = lastIndex+1;
        end

        averageRank = mean(firstIndex:lastIndex);
        ranks(sortOrder(firstIndex:lastIndex)) = averageRank;
        firstIndex = lastIndex+1;
    end
end

function correlation = pearson_base(x,y)

    x = x(:)-mean(x);
    y = y(:)-mean(y);
    denominator = sqrt(sum(x.^2)*sum(y.^2));

    if denominator <= eps
        correlation = NaN;
    else
        correlation = sum(x.*y)/denominator;
        correlation = min(max(correlation,-1),1);
    end
end

function adjustedP = benjamini_hochberg_base(rawP)

    rawP = rawP(:);
    adjustedP = nan(size(rawP));
    valid = isfinite(rawP);
    validP = rawP(valid);

    if isempty(validP)
        return
    end

    [sortedP,sortOrder] = sort(validP);
    numberOfTests = numel(sortedP);
    sortedAdjusted = sortedP.*numberOfTests./(1:numberOfTests)';

    % Enforce monotonicity from the largest rank to the smallest.
    for rankNumber = numberOfTests-1:-1:1
        sortedAdjusted(rankNumber) = min( ...
            sortedAdjusted(rankNumber),sortedAdjusted(rankNumber+1));
    end

    sortedAdjusted = min(sortedAdjusted,1);
    unsortedAdjusted = nan(numberOfTests,1);
    unsortedAdjusted(sortOrder) = sortedAdjusted;
    adjustedP(valid) = unsortedAdjusted;
end

function style_questionnaire_axes(targetAxes)

    box(targetAxes,'off');
    set(targetAxes, ...
        'Color','w', ...
        'FontName','Arial', ...
        'FontSize',10, ...
        'LineWidth',1.4, ...
        'TickDir','out', ...
        'Layer','top');
end

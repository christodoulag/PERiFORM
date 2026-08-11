%% WP4 PPS CLEAN SCRIPT: group plots + slope/R^2
clear; close all; clc;

dbstop if error
% dbstop if warning

%% -------- DATA AND PARAMETERS --------

% Data folder
dataFolder = '';

% Participant identifiers
% The underscore is included because filenames follow: 01AM_pre.csv
name = {''};

% Condition file patterns
strings = {'pre*.csv','post*.csv'};
% User-facing labels. File patterns remain pre/post because they correspond
% to the existing CSV filenames.
condNames = {'Baseline','Post-training'};

% Distances used in the PPS task
distances = [0.25 0.75 1.25 1.75 2.25]';

% Parameters
centralonly = 0;  % 0 = central only, 1 = all, 2 = peripheral only
cutoff = 3;       % RT outlier threshold in standard deviations

% Dimensions
nSubjects = numel(name);
nConditions = numel(strings);
nDistances = numel(distances);

% Mean RT storage
% One row per participant and one column per distance
rt_pre = nan(nSubjects,nDistances);
rt_post = nan(nSubjects,nDistances);

% Linear baseline QC measures calculated for the baseline assessment only
slope_pre = nan(nSubjects,1);
rsq_prelin = nan(nSubjects,1);

% Number of valid trials per participant, condition, and distance
nValidTrials = nan(nSubjects,nConditions,nDistances);

%% -------- LOOP SUBJECTS --------

for s = 1:nSubjects

    % Each cell contains:
    % [Vibrate_distance, RT, Monster]
    sub = cell(1,nConditions);

    %% Load baseline and Post-training files
    for c = 1:nConditions

        files = dir(fullfile( ...
            dataFolder, ...
            [name{s} strings{c}]));

        if isempty(files)

            warning( ...
                'Missing file: %s%s', ...
                name{s},strings{c});

            continue
        end

        if numel(files) > 1

            error( ...
                ['More than one file found for %s%s. ' ...
                 'Keep only one file per condition.'], ...
                 name{s},strings{c});
        end

        filePath = fullfile( ...
            files(1).folder, ...
            files(1).name);

        T = readtable( ...
            filePath, ...
            'VariableNamingRule','preserve');

        % Confirm that the required variables exist
        requiredVariables = ...
            {'Vibrate_distance','RT','Monster'};

        if ~all(ismember( ...
                requiredVariables, ...
                T.Properties.VariableNames))

            error( ...
                ['Required variables are missing from %s. ' ...
                 'Expected: Vibrate_distance, RT, Monster.'], ...
                 files(1).name);
        end

        sub{c} = [ ...
            T.Vibrate_distance, ...
            T.RT, ...
            T.Monster];
    end

    %% Clean each condition
    def = cell(1,nConditions);

    for c = 1:nConditions

        X = sub{c};

        if isempty(X)
            continue
        end

        % Select central or peripheral stimuli
        if centralonly == 0

            X = X(X(:,3) == 0,:);

        elseif centralonly == 2

            X = X(X(:,3) ~= 0,:);
        end

        % Remove catch and tactile-only/control trials
        X(X(:,1) == -1,2) = NaN;
        X(X(:,1) == 0, 2) = NaN;

        % Remove RTs shorter than 100 ms
        % or longer than 1000 ms
        X(X(:,2) < 0.1,2) = NaN;
        X(X(:,2) > 1,  2) = NaN;

        % Remove rows containing invalid RTs
        X = X(~isnan(X(:,2)),:);

        def{c} = X;
    end

    %% Compute mean RT per distance
    subjMeans = nan(nConditions,nDistances);

    for c = 1:nConditions

        X = def{c};

        if isempty(X)
            continue
        end

        for d = 1:nDistances

            trials = X( ...
                X(:,1) == distances(d),2);

            if isempty(trials)
                continue
            end

            % Remove RT outliers within participant,
            % condition, and distance
            mu = mean(trials,'omitnan');
            sd = std(trials,'omitnan');

            if isfinite(sd) && sd > 0

                z = abs((trials-mu)/sd);
                trials = trials(z < cutoff);
            end

            % Store number of valid trials
            nValidTrials(s,c,d) = ...
                sum(isfinite(trials));

            % Store mean RT
            subjMeans(c,d) = ...
                mean(trials,'omitnan');
        end
    end

    %% Store participant-level mean RTs
    rt_pre(s,:) = subjMeans(1,:);
    rt_post(s,:) = subjMeans(2,:);

    %% Calculate linear slope and R² for the baseline assessment only
    % Same baseline QC procedure as WP3
    yPre = rt_pre(s,:)';

    if all(isfinite(yPre)) && ...
            max(yPre) > min(yPre)

        pPre = polyfit(distances,yPre,1);
        yhatPre = polyval(pPre,distances);

        ssresPre = sum((yPre-yhatPre).^2);
        sstotPre = sum((yPre-mean(yPre)).^2);

        if sstotPre > 0

            slope_pre(s) = pPre(1);

            rsq_prelin(s) = ...
                1-(ssresPre/sstotPre);
        end
    end
end

%% -------- IDENTIFY SUBJECTS FAILING BASELINE QC --------

% Criterion 1: Low or invalid baseline R²
lowOrInvalidR2 = ...
    rsq_prelin < 0.10 | ~isfinite(rsq_prelin);

% Criterion 2: Non-positive or invalid baseline slope
nonPositiveSlope = ...
    slope_pre <= 0 | ~isfinite(slope_pre);

% Exclude if either criterion is met
badSubs = lowOrInvalidR2 | nonPositiveSlope;
goodSubs = ~badSubs;

nBad = sum(badSubs);

% Record whether all five Post-training distances are available
postComplete = all(isfinite(rt_post),2);

%% Create an exclusion reason for every participant

exclusionReason = ...
    repmat({'Included'},nSubjects,1);

for i = 1:nSubjects

    if ~isfinite(rsq_prelin(i)) || ...
            ~isfinite(slope_pre(i))

        exclusionReason{i} = ...
            'Invalid baseline slope or R2';

    elseif lowOrInvalidR2(i) && ...
            nonPositiveSlope(i)

        exclusionReason{i} = ...
            'Baseline R2 < 0.10 and non-positive slope';

    elseif lowOrInvalidR2(i)

        exclusionReason{i} = ...
            'Baseline R2 < 0.10';

    elseif nonPositiveSlope(i)

        exclusionReason{i} = ...
            'Non-positive baseline slope';
    end
end

%% Save participant identities before removing anyone

excludedSubjects = name(badSubs);
includedSubjects = name(goodSubs);

%% Create and save QC summary

qcTable = table( ...
    string(name(:)), ...
    slope_pre, ...
    rsq_prelin, ...
    postComplete, ...
    lowOrInvalidR2, ...
    nonPositiveSlope, ...
    badSubs, ...
    exclusionReason, ...
    'VariableNames', ...
    {'Subject', ...
     'BaselineSlope', ...
     'BaselineR2', ...
     'PostTrainingComplete', ...
     'LowOrInvalid_BaselineR2', ...
     'NonPositive_BaselineSlope', ...
     'Excluded_BaselineQC', ...
     'ExclusionReason'});

writetable( ...
    qcTable, ...
    fullfile( ...
        dataFolder, ...
        'WP4_PPS_baseline_QC_summary.csv'));

%% Display QC results

fprintf( ...
    ['\nExcluding %d participant(s) based on ' ...
     'the baseline PPS quality-control criteria:\n'], ...
     nBad);

if nBad > 0

    disp('--------------------------------------------------------------');
    disp(' Subject       Baseline slope       Baseline R²       Reason');
    disp('--------------------------------------------------------------');

    for i = find(badSubs)'

        fprintf( ...
            ' %-10s   %14.4f      %11.3f       %s\n', ...
            name{i}, ...
            slope_pre(i), ...
            rsq_prelin(i), ...
            exclusionReason{i});
    end

    disp('--------------------------------------------------------------');

else

    disp('No participants met the baseline exclusion criteria.');
end

%% -------- FIGURE 1: DETAILED BASELINE QC --------

% Plot all participants before exclusions
figure(1); clf; hold on;

set(gcf,'Color','w');

hIncluded = scatter( ...
    slope_pre(goodSubs), ...
    rsq_prelin(goodSubs), ...
    70,'k','filled');

hExcluded = scatter( ...
    slope_pre(badSubs), ...
    rsq_prelin(badSubs), ...
    70,'r','filled');

% Display both baseline QC thresholds
hR2 = yline( ...
    0.10,'--r', ...
    'R² threshold = 0.10', ...
    'LineWidth',1.5);

hSlope = xline( ...
    0,'--b', ...
    'Slope threshold = 0', ...
    'LineWidth',1.5);

xlabel('Regression slope');
ylabel('Baseline R^2');
title('WP4 baseline PPS quality control');

box off;

set(gca, ...
    'LineWidth',2, ...
    'FontSize',12, ...
    'TickDir','out');

legend( ...
    [hIncluded,hExcluded,hR2,hSlope], ...
    {'Included', ...
     'Excluded', ...
     'R² threshold', ...
     'Slope threshold'}, ...
    'Location','best');

%% -------- FIGURE 2: WP3-STYLE SLOPE vs R² --------

% All participants are displayed before exclusions
figure(2); clf; hold on;

set(gcf, ...
    'Color','w', ...
    'Position',[100 100 640 480]);

validForPlot = ...
    isfinite(slope_pre) & isfinite(rsq_prelin);

% Open circles, consistent with the WP3 figure
plot( ...
    slope_pre(validForPlot), ...
    rsq_prelin(validForPlot), ...
    'o', ...
    'LineStyle','none', ...
    'MarkerSize',7, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',[0.35 0.35 0.35], ...
    'LineWidth',1);

% Display both QC thresholds
yline( ...
    0.10,'--k', ...
    'LineWidth',1.2, ...
    'HandleVisibility','off');

xline( ...
    0,'--k', ...
    'LineWidth',1.2, ...
    'HandleVisibility','off');

xlabel('Regression Slope');
ylabel('R^2');

% Axis limits similar to the WP3 figure
finiteSlopes = slope_pre(validForPlot);

if isempty(finiteSlopes)

    xlim([-0.01 0.05]);

else

    lowerX = min( ...
        -0.01, ...
        min(finiteSlopes)-0.005);

    upperX = max( ...
        0.05, ...
        max(finiteSlopes)+0.005);

    xlim([lowerX upperX]);
end

ylim([0 1]);

box off;

set(gca, ...
    'LineWidth',2, ...
    'FontSize',12, ...
    'TickDir','out');

%% -------- REMOVE SUBJECTS FAILING BASELINE QC --------

rt_pre(badSubs,:) = [];
rt_post(badSubs,:) = [];

nValidTrials(badSubs,:,:) = [];

slope_pre(badSubs) = [];
rsq_prelin(badSubs) = [];

name(badSubs) = [];

% Update participant count after exclusions
nSubjects = numel(name);

%% -------- WP3-STYLE BASELINE AND POST-TRAINING FIGURES --------

% Same y-axis limits used in the WP3 figures
commonYLim = [0.19 0.45];

% Reproducible horizontal jitter
rng(1);

% Figure 3: Baseline assessment
plot_bar_with_dots( ...
    3, ...
    distances, ...
    rt_pre, ...
    'BASELINE', ...
    [0.85 0.20 0.20], ...
    commonYLim, ...
    1.27);

% Figure 4: Post-training assessment following Pull Training
plot_bar_with_dots( ...
    4, ...
    distances, ...
    rt_post, ...
    'POST-TRAINING', ...
    [0.20 0.20 0.85], ...
    commonYLim, ...
    1.68);

%% -------- SAVE CLEANED WORKSPACE --------

% This command must appear before the local function. Save only the cleaned
% numerical and tabular variables so that graphics handles are not written
% into the MAT file.
save( ...
    fullfile( ...
        dataFolder, ...
        'WP4_PPS_clean_workspace.mat'), ...
    'name','nSubjects','distances','rt_pre','rt_post', ...
    'nValidTrials','slope_pre','rsq_prelin','excludedSubjects', ...
    'includedSubjects','qcTable','cutoff','centralonly','dataFolder');

%% -------- HELPER FUNCTION --------
% Local functions must remain at the end of the script

function plot_bar_with_dots( ...
    figNum,distances,rt_mat,titleTxt,barColor,yLimits,ppsBoundary)

    figure(figNum); clf;

    set(gcf, ...
        'Color','w', ...
        'Position',[100 100 640 480]);

    %% Group descriptive statistics
    av = mean(rt_mat,1,'omitnan');
    sd = std(rt_mat,0,1,'omitnan');
    n = sum(isfinite(rt_mat),1);
    sem = sd./sqrt(n);

    valid = isfinite(av) & n > 0;

    %% Group mean bars
    b = bar( ...
        distances(valid), ...
        av(valid), ...
        'LineWidth',1.5);

    b.FaceColor = barColor;
    b.BarWidth = 0.80;

    hold on;
    
    % Display the mean participant-level PPS-boundary estimate.
    xline( ...
        ppsBoundary, ...
        '--k', ...
        sprintf('Mean PPS boundary = %.2f m', ppsBoundary), ...
        'LineWidth', 1.2, ...
        'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'left');
    %% Standard error bars
    errorbar( ...
        distances(valid), ...
        av(valid), ...
        sem(valid), ...
        'k', ...
        'LineWidth',1.5, ...
        'LineStyle','none', ...
        'CapSize',8);

    %% Individual participant values
    xj = distances' + ...
        0.025*randn(size(rt_mat));

    plot( ...
        xj, ...
        rt_mat, ...
        'o', ...
        'LineStyle','none', ...
        'MarkerSize',3.5, ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor',[0.30 0.30 0.30], ...
        'LineWidth',0.6);

    %% Labels and formatting
    title( ...
        titleTxt, ...
        'FontWeight','bold');

    xlabel('Distance (m)');
    ylabel('RTs (s)');

    xlim([0 2.5]);
    ylim(yLimits);

    xticks(distances);
    xticklabels(compose('%.2f',distances));

    yticks(0.20:0.05:0.45);

    box off;

    set(gca, ...
        'LineWidth',1.5, ...
        'FontSize',11, ...
        'TickDir','out');
end


%% -------- FIGURE 2: BASELINE SLOPE vs R² QUALITY CONTROL --------

figure(2);
clf;
hold on;
plot_bar_with_dots
set(gcf, ...
    'Color', 'w', ...
    'Position', [100 100 680 520]);

% Use the QC table because it contains all participants,
% including the participant excluded from subsequent analyses.
slopeAll   = qcTable.BaselineSlope(:);
r2All      = qcTable.BaselineR2(:);
excludedQC = logical(qcTable.Excluded_BaselineQC(:));

% Retain participants with finite values.
validForPlot = isfinite(slopeAll) & isfinite(r2All);

includedForPlot = validForPlot & ~excludedQC;
excludedForPlot = validForPlot & excludedQC;

% Included participants: open grey circles.
hIncluded = scatter( ...
    slopeAll(includedForPlot), ...
    r2All(includedForPlot), ...
    65, ...
    'o', ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [0.35 0.35 0.35], ...
    'LineWidth', 1.2);

% Excluded participant: filled red diamond.
hExcluded = scatter( ...
    slopeAll(excludedForPlot), ...
    r2All(excludedForPlot), ...
    85, ...
    'd', ...
    'MarkerFaceColor', [0.85 0.20 0.20], ...
    'MarkerEdgeColor', [0.60 0.05 0.05], ...
    'LineWidth', 1.3);

% Predefined quality-control thresholds.
yline( ...
    0.10, ...
    '--k', ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

xline( ...
    0, ...
    '--k', ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

% Axis labels.
xlabel('Regression Slope');
ylabel('R^2');

% Axis limits.
finiteSlopes = slopeAll(validForPlot);

if isempty(finiteSlopes)
    xlim([-0.01 0.05]);
else
    lowerX = min(-0.01, min(finiteSlopes) - 0.005);
    upperX = max(0.05, max(finiteSlopes) + 0.005);
    xlim([lowerX upperX]);
end

ylim([0 1]);

% Legend.
legend( ...
    [hIncluded, hExcluded], ...
    {'Included', 'Excluded by baseline QC'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

% Formatting.
box off;

set(gca, ...
    'LineWidth', 1.5, ...
    'FontSize', 12, ...
    'TickDir', 'out');

hold off;

% Export.
exportgraphics( ...
    gcf, ...
    fullfile(dataFolder, ...
    'Figure_6_Baseline_PPS_Quality_Control.png'), ...
    'Resolution', 300);

savefig( ...
    gcf, ...
    fullfile(dataFolder, ...
    'Figure_6_Baseline_PPS_Quality_Control.fig'));
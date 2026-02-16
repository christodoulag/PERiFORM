%% PPS CLEAN SCRIPT (group plots + slope/R^2)
clear; close all; clc;

% IMPORTANT: do NOT stop on warning (you have warnings when files missing)
dbstop if error
% dbstop if warning

%Folder %ADD DIRECTORY
folder=cd('');


% %Subjects %ADD SUBJECTNAMES
name={''};

%Participants excluded
%

% Condition file patterns
strings = {'pre*.csv','tool*.csv','torso*.csv','whole*.csv'};
condNames = {'pre','tool','torso','whole'};

% Distances used in your task
distances = [0.25 0.75 1.25 1.75 2.25]';

% Parameters
centralonly = 0;   % 0 central only, 1 all, 2 peripheral only
cutoff = 3;        % outlier cutoff (z-score)

% Storage: one row per subject, 5 columns (distances)
rt_pre   = nan(length(name),5);
rt_tool  = nan(length(name),5);
rt_torso = nan(length(name),5);
rt_whole = nan(length(name),5);

% For Figure 1 (one dot per subject)
slope_pre  = nan(length(name),1);
rsq_prelin = nan(length(name),1);

%% LOOP SUBJECTS
for s = 1:length(name)

    % Load each condition for this subject
    sub = cell(1,4);  % {pre, tool, torso, whole}, each is Nx3 [distance RT monster]

    for c = 1:4
        files = dir([name{s} strings{c}]);

        if isempty(files)
            warning('Missing: %s%s', name{s}, strings{c});
            continue
        end

        if numel(files) > 1
            error('More than 1 file found for %s%s. Keep only one.', name{s}, strings{c});
        end

        T = readtable(files(1).name, 'VariableNamingRule','preserve');

        sub{c} = [T.Vibrate_distance, T.RT, T.Monster];
    end

    % Clean each condition
    def = cell(1,4);
    for c = 1:4
        X = sub{c};
        if isempty(X), continue; end

        % monster filter
        if centralonly==0
            X = X(X(:,3)==0,:);
        elseif centralonly==2
            X = X(X(:,3)~=0,:);
        end

        % remove catch/base by making RT NaN
        X(X(:,1)==-1,2) = NaN;
        X(X(:,1)==0, 2) = NaN;

        % remove too fast / too slow
        X(X(:,2)<0.1,2) = NaN;
        X(X(:,2)>1,  2) = NaN;

        % drop NaNs
        X = X(~isnan(X(:,2)),:);

        def{c} = X;
    end

    % Compute mean RT per distance with outlier removal
    subjMeans = nan(4,5); % row=condition, col=distance

    for c = 1:4
        X = def{c};
        if isempty(X), continue; end

        for d = 1:5
            trials = X(X(:,1)==distances(d), 2);

            if isempty(trials)
                subjMeans(c,d) = NaN;
                continue
            end

            mu = mean(trials,'omitnan');
            sd = std(trials,'omitnan');

            if sd > 0
                z = abs((trials - mu) / sd);
                trials = trials(z < cutoff);
            end

            subjMeans(c,d) = mean(trials,'omitnan');
        end
    end

    % Save subject means into matrices (one row per subject!)
    rt_pre(s,:)   = subjMeans(1,:);
    rt_tool(s,:)  = subjMeans(2,:);
    rt_torso(s,:) = subjMeans(3,:);
    rt_whole(s,:) = subjMeans(4,:);

    % Figure 1 values (from PRE only): slope + R^2 across the 5 distances
    y = rt_pre(s,:)';
    if all(~isnan(y))
        p = polyfit(distances, y, 1);     % linear fit
        slope_pre(s) = p(1);
        yhat = polyval(p, distances);
        ssres = sum((y - yhat).^2);
        sstot = sum((y - mean(y)).^2);
        rsq_prelin(s) = 1 - ssres/sstot;
    end
end

%% -------- REMOVE SUBJECTS WITH LOW R^2 --------
badSubs = rsq_prelin < 0.1 | isnan(rsq_prelin);  % mark invalid subjects
nBad = sum(badSubs);

fprintf('\nExcluding %d subjects with R^2 < 0.1 or invalid values:\n', nBad);

if nBad > 0
    disp('--------------------------------------------');
    disp(' Subject    R^2 value');
    disp('--------------------------------------------');
    for i = find(badSubs)'
        fprintf(' %-8s   %.3f\n', name{i}, rsq_prelin(i));
    end
    disp('--------------------------------------------');
end

% Remove them from all data matrices
rt_pre(badSubs,:)   = [];
rt_tool(badSubs,:)  = [];
rt_torso(badSubs,:) = [];
rt_whole(badSubs,:) = [];
slope_pre(badSubs)  = [];
rsq_prelin(badSubs) = [];
name(badSubs)       = [];

%% -------- FIGURE 1: slope vs R^2 --------
figure(1);
plot(slope_pre, rsq_prelin, 'ko', 'MarkerSize', 8);
xlabel('Regression Slope'); ylabel('R^2');
ylim([0 1]); box off;
set(gca,'LineWidth',2,'FontSize',12);

%% -------- BAR + DOT PLOTS (Figures 2–5) --------
plot_bar_with_dots(2, distances, rt_pre,   'BASELINE',   [0.85 0.20 0.20]); % red
plot_bar_with_dots(3, distances, rt_tool,  'PULL',  [0.20 0.20 0.85]); % blue
plot_bar_with_dots(4, distances, rt_torso, 'TORSO', [0.20 0.70 0.20]); % green
plot_bar_with_dots(5, distances, rt_whole, 'WHOLE', [0.60 0.20 0.60]); % purple

%% -------- Helper function (at END of script) --------
function plot_bar_with_dots(figNum, distances, rt_mat, titleTxt, barColor)

    figure(figNum); clf;

    av  = mean(rt_mat, 1, 'omitnan');
    sd  = std(rt_mat,  0, 1, 'omitnan');
    n   = sum(~isnan(rt_mat),1);
    sem = sd ./ sqrt(n);

    valid = ~isnan(av);

    b = bar(distances(valid), av(valid), 'LineWidth', 2);
    b.FaceColor = barColor;
    hold on

    errorbar(distances(valid), av(valid), sem(valid), 'k', 'LineWidth', 2, 'LineStyle','none');

    xj = distances' + 0.03*randn(size(rt_mat));
    plot(xj, rt_mat, 'ko', 'MarkerSize', 5);

    title(titleTxt);
    xlabel('Distance (m)'); ylabel('RTs (s)');

    xlim([0 2.5]);
    ylim([0.19 0.45]);

    % add this
    xticks(distances);
    xticklabels(compose('%.2f', distances));

    box off;
    set(gca,'LineWidth',2,'FontSize',12);
end


save('C:\Users\chris\Desktop\workspace.mat')
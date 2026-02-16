%% Gender secondary analysis (ONE FILE, FIXED + DETAILED BOOTSTRAP PROGRESS)
% Shows EXACTLY where you are during bootstrapping by printing:
% - attempts counter (every cfg.attemptEvery attempts)
% - accepted counter (every cfg.progressEvery accepted)
% - last R^2 and last x0 (PPS border) from the most recent fit attempt
% - elapsed time since the start of each condition
%
% IMPORTANT:
% - original code only printed when ACCEPTED samples reached 250.
%   If acceptance is slow (R^2<0.8 often), you see "nothing".
% - This version prints progress based on ATTEMPTS too.
%
% NOTE: Replaced BOXplot (needs Stats toolbox) with a toolbox-free plot.

%% Clean the data
clear all;
close all;
clc;

warning('on','all');

%% Folder %ADD DIRECTORY
folderPath = '';
cd(folderPath);

%% Subject lists %ADD SUBJECTNAMES
women_names = {''};
men_names   = {''};

%Participants excluded
%
%% Parameters
cfg = struct();
cfg.centralonly   = 0;
cfg.cutoff        = 3;
cfg.fig_per_group = 0;
cfg.howmanyboot   = 12001;
cfg.waitfit       = 0;
cfg.progressEvery = 250;    % accepted prints
cfg.attemptEvery  = 250;    % attempt prints (always shows progress)

%% Run WOMEN then MEN
fprintf('\nStarting WOMEN analysis...\n');
tic;
resW = run_group_pps(women_names, folderPath, cfg, 'WOMEN');
fprintf('Finished WOMEN in %.1f seconds.\n', toc);

fprintf('\nStarting MEN analysis...\n');
tic;
resM = run_group_pps(men_names, folderPath, cfg, 'MEN');
fprintf('Finished MEN in %.1f seconds.\n', toc);

%% Compare (Women - Men)
fprintf('\n================ GENDER COMPARISON (Women - Men) ================\n');
fprintf('N women=%d | N men=%d\n', resW.Nsubs, resM.Nsubs);

report_diff('PRE',   resW.BootBorderSigPre,   resM.BootBorderSigPre);
report_diff('TOOL',  resW.BootBorderSigTool,  resM.BootBorderSigTool);
report_diff('TORSO', resW.BootBorderSigTorso, resM.BootBorderSigTorso);
report_diff('WHOLE', resW.BootBorderSigWhole, resM.BootBorderSigWhole);
fprintf('=================================================================\n');
%%%% ONE FIGURE: PPS border by gender (toolbox-free "boxplot-like" summary)
conds = {'BASELINE','PULL','TORSO','WHOLE'};

W = {resW.BootBorderSigPre, ...
     resW.BootBorderSigTool, ...
     resW.BootBorderSigTorso, ...
     resW.BootBorderSigWhole};

M = {resM.BootBorderSigPre, ...
     resM.BootBorderSigTool, ...
     resM.BootBorderSigTorso, ...
     resM.BootBorderSigWhole};

figure('Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

%% ===== GLOBAL Y-LIMITS (same for all subplots) =====
allDataGlobal = [ ...
    W{1}(:); W{2}(:); W{3}(:); W{4}(:); ...
    M{1}(:); M{2}(:); M{3}(:); M{4}(:)];

allDataGlobal = allDataGlobal(~isnan(allDataGlobal));

globalMin = min(allDataGlobal);
globalMax = max(allDataGlobal);

padding = 0.05 * (globalMax - globalMin);

ymin = globalMin - padding;
ymax = globalMax + padding;

% optional rounding for nicer axis
ymin = floor(ymin*10)/10;
ymax = ceil(ymax*10)/10;

%% ===== PLOT LOOP =====
for i = 1:4
    nexttile; hold on;

    % jittered scatter
    jitterW = (rand(size(W{i})) - 0.5) * 0.15;
    jitterM = (rand(size(M{i})) - 0.5) * 0.15;

    xW = 1 + jitterW;
    xM = 2 + jitterM;

    plot(xW, W{i}, '.', 'MarkerSize', 6);
    plot(xM, M{i}, '.', 'MarkerSize', 6);

    % summary stats
    draw_summary(1, W{i});
    draw_summary(2, M{i});

    % axis settings
    xlim([0.5 2.5]);
    ylim([ymin ymax]);

    set(gca,'XTick',[1 2],'XTickLabel',{'Women','Men'});

    title(conds{i});
    ylabel('PPS border (m)');

    set(gca,'FontSize',14,'LineWidth',1.5);
    box off;
end

sgtitle('PPS Border by Gender (bootstrap distributions)', ...
        'FontSize',16,'FontWeight','bold');

% %% ONE FIGURE: PPS border by gender (toolbox-free "boxplot-like" summary)
% conds = {'BASELINE','PULL','TORSO','WHOLE'};
% 
% W = {resW.BootBorderSigPre, ...
%      resW.BootBorderSigTool, ...
%      resW.BootBorderSigTorso, ...
%      resW.BootBorderSigWhole};
% 
% M = {resM.BootBorderSigPre, ...
%      resM.BootBorderSigTool, ...
%      resM.BootBorderSigTorso, ...
%      resM.BootBorderSigWhole};
% 
% figure('Color','w');
% tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
% 
% for i = 1:4
%     nexttile; hold on;
% 
%     % jittered scatter (distributions)
%     jitterW = (rand(size(W{i})) - 0.5) * 0.15;
%     jitterM = (rand(size(M{i})) - 0.5) * 0.15;
% 
%     xW = 1 + jitterW;
%     xM = 2 + jitterM;
% 
%     plot(xW, W{i}, '.', 'MarkerSize', 6);
%     plot(xM, M{i}, '.', 'MarkerSize', 6);
% 
%     % summary stats (median + IQR whiskers)
%     draw_summary(1, W{i});
%     draw_summary(2, M{i});
% 
%     ax = gca;
%     ax.XLim = [0 2.5];
%     set(gca,'XTick',[1 2],'XTickLabel',{'Women','Men'});
% 
%     title(conds{i});
%     ylabel('PPS border (m)');
%     %ylim([0.8 2.2]);
%     allData = [W{i}; M{i}];
% allData = allData(~isnan(allData));
% 
% if ~isempty(allData)
%     ymin = min(allData);
%     ymax = max(allData);
% 
%     if ymax > ymin
%         padding = 0.05 * (ymax - ymin);
%         ylim([ymin - padding, ymax + padding]);
%     end
% end
% 
%     set(gca,'FontSize',14,'LineWidth',1.5);
%     box off;
% end
% 
% sgtitle('PPS Border by Gender (bootstrap distributions)', ...
%         'FontSize',16,'FontWeight','bold');

%% ========================= LOCAL FUNCTIONS =========================

function out = run_group_pps(name, folderPath, cfg, suffix)

    cd(folderPath);

    eqsig = '1/(1+exp(-(x-x0)/b))';

    centralonly   = cfg.centralonly;
    cutoff        = cfg.cutoff;
    howmanyboot   = cfg.howmanyboot;
    progressEvery = cfg.progressEvery;
    attemptEvery  = cfg.attemptEvery;

    strings = {'pre*.csv','tool*.csv','torso*.csv', 'whole*.csv'};

    supsubdatapre   = [];
    supsubdatatool  = [];
    supsubdatatorso = [];
    supsubdatawhole = [];

    fprintf('[%s] Reading files...\n', suffix);
    for l = 1:length(name)
        for g = 1:length(strings)
            files = dir(strcat(name{l}, strings{g}));
            if isempty(files)
                warning('[%s] Missing file(s) for %s %s', suffix, name{l}, strings{g});
                continue
            end
            data = readtable(files(1).name, 'VariableNamingRule', 'preserve');

            if g==1
                supsubdatapre   = [supsubdatapre;   data.Vibrate_distance, data.RT, data.Monster]; %#ok<AGROW>
            elseif g==2
                supsubdatatool  = [supsubdatatool;  data.Vibrate_distance, data.RT, data.Monster]; %#ok<AGROW>
            elseif g==3
                supsubdatatorso = [supsubdatatorso; data.Vibrate_distance, data.RT, data.Monster]; %#ok<AGROW>
            elseif g==4
                supsubdatawhole = [supsubdatawhole; data.Vibrate_distance, data.RT, data.Monster]; %#ok<AGROW>
            end
        end
    end

    if isempty(supsubdatapre) || isempty(supsubdatatool) || isempty(supsubdatatorso) || isempty(supsubdatawhole)
        error('[%s] One or more conditions have NO DATA after file reading. Check filenames/folderPath.', suffix);
    end

    ppsdata = {supsubdatapre, supsubdatatool, supsubdatatorso, supsubdatawhole};
    dist    = [0.25; 0.75; 1.25; 1.75; 2.25];

    fprintf('[%s] Cleaning data...\n', suffix);
    for k = 1:length(strings)

        if centralonly==0
            central_data = find(ppsdata{k}(:,3)==0);
            ppsdata{k}   = ppsdata{k}(central_data,:);
        elseif centralonly==2
            central_data = find(ppsdata{k}(:,3)~=0);
            ppsdata{k}   = ppsdata{k}(central_data,:);
        end

        catcht = find(ppsdata{k}(:,1)==-1);
        ppsdata{k}(catcht,2:3) = NaN;
        base   = find(ppsdata{k}(:,1)==0);
        ppsdata{k}(base,2:3)   = NaN;
        missed = find(ppsdata{k}(:,2)<0.1);
        ppsdata{k}(missed,:)   = NaN;
        longt  = find(ppsdata{k}(:,2)>1);
        ppsdata{k}(longt,:)    = NaN;

        A = isnan(ppsdata{k}(:,2));
        ppsdata{k} = ppsdata{k}(~A,:);

        distance1 = find(ppsdata{k}(:,1)==0.25);
        distance2 = find(ppsdata{k}(:,1)==0.75);
        distance3 = find(ppsdata{k}(:,1)==1.25);
        distance4 = find(ppsdata{k}(:,1)==1.75);
        distance5 = find(ppsdata{k}(:,1)==2.25);

        if k==1
            predist1=ppsdata{k}(distance1,:);
            predist2=ppsdata{k}(distance2,:);
            predist3=ppsdata{k}(distance3,:);
            predist4=ppsdata{k}(distance4,:);
            predist5=ppsdata{k}(distance5,:);
        elseif k==2
            tooldist1=ppsdata{k}(distance1,:);
            tooldist2=ppsdata{k}(distance2,:);
            tooldist3=ppsdata{k}(distance3,:);
            tooldist4=ppsdata{k}(distance4,:);
            tooldist5=ppsdata{k}(distance5,:);
        elseif k==3
            torsodist1=ppsdata{k}(distance1,:);
            torsodist2=ppsdata{k}(distance2,:);
            torsodist3=ppsdata{k}(distance3,:);
            torsodist4=ppsdata{k}(distance4,:);
            torsodist5=ppsdata{k}(distance5,:);
        elseif k==4
            wholedist1=ppsdata{k}(distance1,:);
            wholedist2=ppsdata{k}(distance2,:);
            wholedist3=ppsdata{k}(distance3,:);
            wholedist4=ppsdata{k}(distance4,:);
            wholedist5=ppsdata{k}(distance5,:);
        end
    end

    optCommon = @( ) fitoptions(eqsig, 'Lower',[eps 0.25], 'Upper',[Inf 2.25], 'StartPoint',[0.5 1.25]);

    [pps_border_presig, ~]   = bootstrap_condition('PRE',   predist1,predist2,predist3,predist4,predist5, dist, eqsig, optCommon, cutoff, howmanyboot, progressEvery, attemptEvery, suffix);
    [pps_border_toolsig, ~]  = bootstrap_condition('TOOL',  tooldist1,tooldist2,tooldist3,tooldist4,tooldist5, dist, eqsig, optCommon, cutoff, howmanyboot, progressEvery, attemptEvery, suffix);
    [pps_border_torsosig, ~] = bootstrap_condition('TORSO', torsodist1,torsodist2,torsodist3,torsodist4,torsodist5, dist, eqsig, optCommon, cutoff, howmanyboot, progressEvery, attemptEvery, suffix);
    [pps_border_wholesig, ~] = bootstrap_condition('WHOLE', wholedist1,wholedist2,wholedist3,wholedist4,wholedist5, dist, eqsig, optCommon, cutoff, howmanyboot, progressEvery, attemptEvery, suffix);

    BootBorderSigPre    = [pps_border_presig{:}];
    BootBorderSigTool   = [pps_border_toolsig{:}];
    BootBorderSigTorso  = [pps_border_torsosig{:}];
    BootBorderSigWhole  = [pps_border_wholesig{:}];

    save(['BootBorderSigPre_' suffix '.mat'],   'BootBorderSigPre');
    save(['BootBorderSigTool_' suffix '.mat'],  'BootBorderSigTool');
    save(['BootBorderSigTorso_' suffix '.mat'], 'BootBorderSigTorso');
    save(['BootBorderSigWhole_' suffix '.mat'], 'BootBorderSigWhole');

    out = struct();
    out.suffix = suffix;
    out.Nsubs  = numel(name);
    out.BootBorderSigPre   = BootBorderSigPre;
    out.BootBorderSigTool  = BootBorderSigTool;
    out.BootBorderSigTorso = BootBorderSigTorso;
    out.BootBorderSigWhole = BootBorderSigWhole;
end

function [pps_border_sig, rsq_sig] = bootstrap_condition(label, d1,d2,d3,d4,d5, dist, eqsig, optCommon, cutoff, howmanyboot, progressEvery, attemptEvery, suffix)

    fprintf('[%s] Bootstrapping %s...\n', suffix, label);
    t0 = tic;

    pps_border_sig = cell(1, howmanyboot-1);
    rsq_sig        = cell(1, howmanyboot-1);

    k = 1;
    attempt = 0;
    last_rsq = NaN;
    last_x0  = NaN;

    while k < howmanyboot
        attempt = attempt + 1;

        if mod(attempt, attemptEvery) == 0
            fprintf('[%s] %s attempts=%d | accepted=%d/%d | last R2=%.3f | last x0=%.3f | elapsed=%.1fs\n', ...
                suffix, label, attempt, k-1, howmanyboot-1, last_rsq, last_x0, toc(t0));
        end

        rng(attempt);

        s1 = d1(randi(size(d1,1), [size(d1,1),1]), :);
        s2 = d2(randi(size(d2,1), [size(d2,1),1]), :);
        s3 = d3(randi(size(d3,1), [size(d3,1),1]), :);
        s4 = d4(randi(size(d4,1), [size(d4,1),1]), :);
        s5 = d5(randi(size(d5,1), [size(d5,1),1]), :);

        keep1 = abs((s1(:,2)-mean(s1(:,2)))/std(s1(:,2))) < cutoff;
        keep2 = abs((s2(:,2)-mean(s2(:,2)))/std(s2(:,2))) < cutoff;
        keep3 = abs((s3(:,2)-mean(s3(:,2)))/std(s3(:,2))) < cutoff;
        keep4 = abs((s4(:,2)-mean(s4(:,2)))/std(s4(:,2))) < cutoff;
        keep5 = abs((s5(:,2)-mean(s5(:,2)))/std(s5(:,2))) < cutoff;

        av = [mean(s1(keep1,2)); mean(s2(keep2,2)); mean(s3(keep3,2)); mean(s4(keep4,2)); mean(s5(keep5,2))];

        den = max(av) - min(av);
        if den == 0
            last_rsq = NaN; last_x0 = NaN;
            continue
        end

        y = (av - min(av)) ./ den;

        try
            [fitobj, goodness] = fit(dist, y, eqsig, optCommon());
            last_rsq = goodness.rsquare;
            last_x0  = fitobj.x0;
        catch
            last_rsq = NaN;
            last_x0  = NaN;
        end

        if ~isnan(last_rsq) && last_rsq >= 0.8
            rsq_sig{k}        = last_rsq;
            pps_border_sig{k} = last_x0;

            if mod(k, progressEvery) == 0
                fprintf('[%s] %s ACCEPTED=%d/%d | attempt=%d | R2=%.3f | x0=%.3f | elapsed=%.1fs\n', ...
                    suffix, label, k, howmanyboot-1, attempt, last_rsq, last_x0, toc(t0));
            end

            k = k + 1;
        end
    end

    fprintf('[%s] Finished %s: accepted=%d/%d in %.1fs (attempts=%d)\n', ...
        suffix, label, howmanyboot-1, howmanyboot-1, toc(t0), attempt);
end

function report_diff(label, bootW, bootM)
    d  = bootW - bootM;
    mu = mean(d);
    ci = prctile(d, [2.5 97.5]);
    p  = 2 * min(mean(d >= 0), mean(d <= 0));
    fprintf('%-5s Δ(W-M)=%.4f | 95%% CI [%.4f %.4f] | p=%.4f\n', label, mu, ci(1), ci(2), p);
end

function draw_summary(x, data)
% toolbox-free "boxplot-like" summary: median line + IQR + whiskers (1.5*IQR)
    data = data(:);
    data = data(~isnan(data));

    if isempty(data)
        return
    end

    % quartiles without toolbox
    s = sort(data);
    q1 = quantile_simple(s, 0.25);
    q2 = quantile_simple(s, 0.50);
    q3 = quantile_simple(s, 0.75);
    iqr = q3 - q1;

    low  = max(min(s), q1 - 1.5*iqr);
    high = min(max(s), q3 + 1.5*iqr);

    % IQR "box"
    plot([x-0.12 x+0.12],[q1 q1],'k-','LineWidth',2);
    plot([x-0.12 x+0.12],[q3 q3],'k-','LineWidth',2);
    plot([x-0.12 x-0.12],[q1 q3],'k-','LineWidth',2);
    plot([x+0.12 x+0.12],[q1 q3],'k-','LineWidth',2);

    % median
    plot([x-0.12 x+0.12],[q2 q2],'k-','LineWidth',3);

    % whiskers
    plot([x x],[low q1],'k-','LineWidth',1.5);
    plot([x x],[q3 high],'k-','LineWidth',1.5);
    plot([x-0.06 x+0.06],[low low],'k-','LineWidth',1.5);
    plot([x-0.06 x+0.06],[high high],'k-','LineWidth',1.5);
end

function q = quantile_simple(sortedData, p)
% simple linear-interpolation quantile for already-sorted vector
    n = numel(sortedData);
    if n == 1
        q = sortedData(1);
        return
    end
    pos = 1 + (n-1)*p;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = sortedData(lo);
    else
        w = pos - lo;
        q = (1-w)*sortedData(lo) + w*sortedData(hi);
    end
end

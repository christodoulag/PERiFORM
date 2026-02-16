clear; clc;
%Folder 
folder=cd('');

%ADD DIRECTORY
folderPath = '';
cd(folderPath);

%ADD SUBJECTNAMES
subjects = [{}];

subjects = subjects(:);

conditions = {'pre','tool','torso','whole'};
dist = [0.25 0.75 1.25 1.75 2.25]';

eqsig = '1/(1+exp(-(x-x0)/b))';
opt = fitoptions(eqsig,...
    'Lower',[eps 0.25],...
    'Upper',[Inf 2.25],...
    'StartPoint',[0.5 1.25]);

results = [];

%% ================= SUBJECT-LEVEL PPS =================

for s = 1:length(subjects)

    subj = subjects{s};
    PPS_values = nan(1,4);

    for c = 1:4

        file = dir([subj conditions{c} '*.csv']);
        if isempty(file)
            continue
        end

        data = readtable(file(1).name);

        % Keep central trials
        data = data(data.Monster==0,:);

        % Clean RT
        data.RT(data.RT<0.1 | data.RT>1) = NaN;
        data = data(~isnan(data.RT),:);

        av = nan(5,1);

        for d = 1:5
            av(d) = mean(data.RT(data.Vibrate_distance==dist(d)));
        end

        if any(isnan(av)) || max(av)-min(av)==0
            PPS_values(c) = NaN;
            continue
        end

        % Normalize
        y = (av - min(av)) ./ (max(av)-min(av));

        try
            fitobj = fit(dist, y, eqsig, opt);
            PPS_values(c) = fitobj.x0;
        catch
            PPS_values(c) = NaN;
        end
    end

    results = [results;
        {subj, PPS_values(1), PPS_values(2), PPS_values(3), PPS_values(4)}];
end

PPS = cell2table(results,...
    'VariableNames',{'Subject','Baseline','Pull','Torso','Whole'});

writetable(PPS,'PPS_subject_level.xlsx');

%% ================= DELTA FROM BASELINE =================

PPS.DeltaPull  = ((PPS.Pull  - PPS.Baseline) ./ PPS.Baseline) * 100;
PPS.DeltaTorso = ((PPS.Torso - PPS.Baseline) ./ PPS.Baseline) * 100;
PPS.DeltaWhole = ((PPS.Whole - PPS.Baseline) ./ PPS.Baseline) * 100;

writetable(PPS,'PPS_with_deltas.xlsx');

%% ================= MERGE WITH QUESTIONNAIRES =================

Q = readtable('Questionnaires.xlsx');

merged = innerjoin(PPS, Q, 'Keys','Subject');

writetable(merged,'Merged_PPS_Questionnaires.xlsx');

data = merged;

%% ================= REGRESSION ANALYSIS =================

questionnaires = {'ASRS','ATTC','SSQ'};
outcomes = {'Baseline','DeltaPull','DeltaTorso','DeltaWhole'};

all_p = [];
results_reg = [];

for q = 1:length(questionnaires)

    for o = 1:length(outcomes)

        X = data.(questionnaires{q});
        Y = data.(outcomes{o});

        valid = ~isnan(X) & ~isnan(Y);
        X = X(valid);
        Y = Y(valid);

        if length(X) < 5
            continue
        end

        % Standardized beta
        Xz = zscore(X);
        Yz = zscore(Y);

        mdl = fitlm(Xz, Yz);
        beta = mdl.Coefficients.Estimate(2);
        pval = mdl.Coefficients.pValue(2);

        % Robust regression
        mdl_rob = fitlm(Xz, Yz, 'RobustOpts','on');
        p_rob = mdl_rob.Coefficients.pValue(2);

        all_p = [all_p; pval];

        results_reg = [results_reg;
            {questionnaires{q}, outcomes{o}, beta, pval, p_rob}];

        % Plot
        figure;
        scatter(X, Y, 60, 'filled');
        lsline;
        xlabel(questionnaires{q});
        ylabel(outcomes{o});
        title([questionnaires{q} ' predicting ' outcomes{o}]);
        set(gca,'FontSize',14,'LineWidth',1.5);
        box off;

    end
end

%% ================= FDR (Benjamini-Hochberg) =================

[p_sorted, sort_idx] = sort(all_p);
m = length(all_p);
adj_p = zeros(size(all_p));

for i = 1:m
    adj_p(sort_idx(i)) = p_sorted(i) * m / i;
end

adj_p(adj_p>1) = 1;

results_table = cell2table(results_reg,...
    'VariableNames',{'Questionnaire','Outcome',...
    'Std_Beta','p_linear','p_robust'});

results_table.FDR_p = adj_p;

writetable(results_table,'Questionnaire_Regression_Results.xlsx');

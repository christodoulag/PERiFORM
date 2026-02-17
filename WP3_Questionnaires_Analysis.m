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

results = {};

%% ================= SUBJECT-LEVEL PPS =================

for s = 1:length(subjects)

    subj = subjects{s};
    subj_clean = erase(subj,'_');   % remove underscore

    PPS_values = nan(1,4);

    for c = 1:4

        file = dir([subj conditions{c} '*.csv']);
        if isempty(file)
            continue
        end

        data = readtable(file(1).name,'VariableNamingRule','preserve');

        data = data(data.Monster==0,:);
        data.RT(data.RT<0.1 | data.RT>1) = NaN;
        data = data(~isnan(data.RT),:);

        av = nan(5,1);

        for d = 1:5
            av(d) = mean(data.RT(data.Vibrate_distance==dist(d)));
        end

        if any(isnan(av)) || max(av)==min(av)
            continue
        end

        y = (av - min(av)) ./ (max(av)-min(av));

        try
            fitobj = fit(dist, y, eqsig, opt);
            PPS_values(c) = fitobj.x0;
        catch
            PPS_values(c) = NaN;
        end
    end

    results(end+1,:) = ...
        {subj_clean, PPS_values(1), PPS_values(2), PPS_values(3), PPS_values(4)};
end

PPS = cell2table(results,...
    'VariableNames',{'Subject','Baseline','Pull','Torso','Whole'});

%% ================= DELTAS =================

PPS.DeltaPull  = ((PPS.Pull  - PPS.Baseline) ./ PPS.Baseline) * 100;
PPS.DeltaTorso = ((PPS.Torso - PPS.Baseline) ./ PPS.Baseline) * 100;
PPS.DeltaWhole = ((PPS.Whole - PPS.Baseline) ./ PPS.Baseline) * 100;

%% ================= MERGE =================

Q = readtable('Questionnaires.xlsx','VariableNamingRule','preserve');

disp('Excel column names detected:')
disp(Q.Properties.VariableNames')

merged = innerjoin(PPS, Q, 'Keys','Subject');

disp('Merged table size:')
disp(size(merged))

data = merged;

%% ================= REGRESSION (NO TOOLBOX) =================

questionnaire_names = {'ASRS','ATTC','SSQ','TEQ'};

disp('Using questionnaire columns:')
disp(questionnaire_names)

outcomes = {'Baseline','DeltaPull','DeltaTorso','DeltaWhole'};

results_reg = {};
all_p = [];

nPerm = 5000;

for q = 1:length(questionnaire_names)

    if ~ismember(questionnaire_names{q}, data.Properties.VariableNames)
        continue
    end

    for o = 1:length(outcomes)

        X = data.(questionnaire_names{q});
        Y = data.(outcomes{o});

        valid = ~isnan(X) & ~isnan(Y);
        X = X(valid);
        Y = Y(valid);

        n = length(X);

        if n < 6
            continue
        end

        % Manual standardization
        Xz = (X - mean(X)) ./ std(X);
        Yz = (Y - mean(Y)) ./ std(Y);

        % OLS regression
        Xmat = [ones(n,1) Xz];
        b = Xmat \ Yz;
        beta = b(2);

        % Permutation test
        permBetas = zeros(nPerm,1);

        for p = 1:nPerm
            Yperm = Yz(randperm(n));
            bperm = Xmat \ Yperm;
            permBetas(p) = bperm(2);
        end

        pval = mean(abs(permBetas) >= abs(beta));

        all_p(end+1,1) = pval;

        results_reg(end+1,:) = ...
            {questionnaire_names{q}, outcomes{o}, beta, pval};

        % Plot
        figure;
        scatter(X,Y,70,'filled'); hold on
        xx = linspace(min(X),max(X),100);
        yy = mean(Y) + beta*((xx-mean(X))/std(X))*std(Y);
        plot(xx,yy,'LineWidth',2)

        xlabel(questionnaire_names{q});
        ylabel(outcomes{o});
        title([questionnaire_names{q} ' predicting ' outcomes{o} ...
               ' | beta=' num2str(beta,2) ...
               ' | p=' num2str(pval,2)]);

        set(gca,'FontSize',14,'LineWidth',1.5);
        box off;
    end
end

%% ================= FDR =================

if isempty(results_reg)
    error('No regression results computed.');
end

[p_sorted, sort_idx] = sort(all_p);
m = length(all_p);
adj_p = zeros(m,1);

for i = 1:m
    adj_p(sort_idx(i)) = p_sorted(i) * m / i;
end

adj_p(adj_p>1) = 1;

results_table = cell2table(results_reg,...
    'VariableNames',{'Questionnaire','Outcome',...
    'Std_Beta','p_perm'});

results_table.FDR_p = adj_p;

writetable(results_table,'Questionnaire_Regression_Results.xlsx');

disp('Analysis complete.');

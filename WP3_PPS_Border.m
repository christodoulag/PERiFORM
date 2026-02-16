% Clean everything
clear all; 
close all;
clc;

%Folder %ADD DIRECTORY
folder=cd('');


% %Subjects %ADD SUBJECTNAMES
name={''};

%Participants excluded
%


%Parameters to change
centralonly=0; %0 only stimuli from the central monster 1 all data 2 periph only

cutoff=3; %cut off std value for outliers

%Parameters for the fit 
eq='(ymin+(ymax*exp((x-x0)/b)))/(1+exp((x-x0)/b))'; %sigmoid equation for the fit taken from Serino et al. 2018

options=fitoptions(eq, 'Lower', [0 0.25 1 0], 'Upper', [Inf 2.25 1 0], 'StartPoint', [Inf 1.25 0.5 0.5]);

%start point order b x0 ymax ymin



%strings that mark each condition
strings={'pre*.csv','tool*.csv','torso*.csv', 'whole*.csv'};


% Take the data from the excel, we only want the distance and rts column and which monster is shooting
%Prepare empty matrix to put them 
supsubdatapre=[];
supsubdatatool=[];
supsubdatatorso=[];
supsubdatawhole=[];

 %data pooled from all subjects and put them in separate matrices
for l=1:length(name)
    
    for g=1:length(strings)
    
    files = dir(strcat(name{l}, strings{g}));

% if no file found, skip
if isempty(files)
    warning('No file found for %s%s', name{l}, strings{g});
    continue
end

% if more than 1 file found, stop (so you see the problem)
if numel(files) > 1
    error('More than 1 file found for %s%s. Please keep only one file per condition.', name{l}, strings{g});
end

% read the single file with full path
data = readtable(fullfile(files(1).folder, files(1).name), ...
                 'VariableNamingRule','preserve');
    if g==1 %pre
    supsubdatapre=[supsubdatapre; data.Vibrate_distance, data.RT, data.Monster];
    elseif g==2 %tool
    supsubdatatool=[supsubdatatool; data.Vibrate_distance, data.RT, data.Monster];
    elseif g==3 %torso
    supsubdatatorso=[supsubdatatorso; data.Vibrate_distance, data.RT, data.Monster];
    elseif g==4 %whole
     supsubdatawhole=[supsubdatawhole; data.Vibrate_distance, data.RT, data.Monster];
 
    end
    
    end
    
end


ppsdata_raw={supsubdatapre, supsubdatatool, supsubdatatorso, supsubdatawhole};

fprintf('\nLoaded rows:\n');
fprintf('pre   : %d\n', size(supsubdatapre,1));
fprintf('tool  : %d\n', size(supsubdatatool,1));
fprintf('torso : %d\n', size(supsubdatatorso,1));
fprintf('whole : %d\n\n', size(supsubdatawhole,1));

%clean raw data 
for i = 1:numel(ppsdata_raw)

    if isempty(ppsdata_raw{i})
        warning('ppsdata_raw{%d} is empty (no data loaded)', i);
        ppsdata_def{i} = [];
        continue
    end
    
%Take only trials from the monster selected at line 21    
if centralonly==0
central_data=find(ppsdata_raw{i}(:,3)==0);
ppsdata_raw{i}=ppsdata_raw{i}(central_data,:);
elseif centralonly==2
central_data=find(ppsdata_raw{i}(:,3)~=0);
ppsdata_raw{i}=ppsdata_raw{i}(central_data,:);
end

catcht{i}=find(ppsdata_raw{i}(:,1)==-1); %take out catch trials 
ppsdata_raw{i}(catcht{i},2:3)=NaN; %Put a NaN there 
base{i}=find(ppsdata_raw{i}(:,1)==0); %take out base trials 
ppsdata_raw{i}(base{i},2:3)=NaN; %Put a NaN there 
missed{i}=find(ppsdata_raw{i}(:,2)<0.1); %take out any missed trial (Rt= -999) and RTs too short 
ppsdata_raw{i}(missed{i},2:3)=NaN; %Put a NaN there 
long{i}=find(ppsdata_raw{i}(:,2)>1); %take out trials longer than 1
ppsdata_raw{i}(long{i},2:3)=NaN; %Put a NaN there
 A{i}=isnan(ppsdata_raw{i}(:,2));
ppsdata_def{i}=ppsdata_raw{i}(~A{i},:); %take out NaNs 
end

distances = [0.25 0.75 1.25 1.75 2.25]'; %how many distances should be 6 (0 0.25 0.75 1.25 1.75 2.25) but NaN are counted once for each so is much longer 

% Bootper pre


    scramble_data_pre=ppsdata_def{1};

    
    %now we finally calculate the means 
    for g=1:5 %for the 5 distances 
    eachdist=find(scramble_data_pre(:,1)==distances(g)); %find trials for that distance
    mytrials=scramble_data_pre(eachdist,2);
    total_dist{g}=mytrials;
    av_rt_pre(g)=mean(mytrials,'omitnan');
    std_rt_pre(g)=std(mytrials,'omitnan');
    zeta=[]; 
    zeta=abs((mytrials-av_rt_pre(g))/std_rt_pre(g));
    mytrials(:,3)=zeta>=cutoff; %writes 1 on the trials to discard 
    trialskeep=find(mytrials(:,3)==0);
    mytrials=mytrials(trialskeep,1);
    howmanyout_pre(g,:)=sum(zeta>cutoff);
    av_rt_pre(g)=mean(mytrials,'omitnan'); %new average for each distance of the cleaned data 
    std_rt_pre(g)=std(mytrials,'omitnan');
    end
    
    
    %we normalize the new averages between 0 and 1
   for l=1:length(av_rt_pre)
   norm_av_pre(l)=(av_rt_pre(l)-min(av_rt_pre))/(max(av_rt_pre)-min(av_rt_pre));
   end
    
  % norm_av_pre=av_rt_pre;

   optpre=fitoptions(eq, 'Lower', [0 0.25 max(norm_av_pre) min(norm_av_pre)], 'Upper', [Inf 2.25 max(norm_av_pre) min(norm_av_pre)], 'StartPoint', [Inf 1.25 max(norm_av_pre)/2  max(norm_av_pre)/2]);

   
try 
    
[fitobj_pre, goodness_pre, output_pre, convmsg_pre]=fit(distances,norm_av_pre',eq,optpre);
   rsq_pre=goodness_pre.rsquare;
   pps_border_pre=fitobj_pre.x0;
   
catch 

   rsq_pre=NaN;
   pps_border_pre=NaN;

end



% Bootper tool


    scramble_data_tool=ppsdata_def{2};

    
    %now we finally calculate the means 
    for g=1:5 %for the 5 distances 
    eachdist=find(scramble_data_tool(:,1)==distances(g)); %find trials for that distance
    mytrials=scramble_data_tool(eachdist,2);
    total_dist{g}=mytrials;
    av_rt_tool(g)=mean(mytrials,'omitnan');
    std_rt_tool(g)=std(mytrials,'omitnan');
    zeta=[]; 
    zeta=abs((mytrials-av_rt_tool(g))/std_rt_tool(g));
    mytrials(:,3)=zeta>=cutoff; %writes 1 on the trials to discard 
    trialskeep=find(mytrials(:,3)==0);
    mytrials=mytrials(trialskeep,1);
    howmanyout_tool(g,:)=sum(zeta>cutoff);
    av_rt_tool(g)=mean(mytrials,'omitnan');%new average for each distance of the cleaned data 
    std_rt_tool(g)=std(mytrials,'omitnan');
    end
    
    
    %we normalize the new averages between 0 and 1
   for l=1:length(av_rt_tool)
   norm_av_tool(l)=(av_rt_tool(l)-min(av_rt_tool))/(max(av_rt_tool)-min(av_rt_tool));
   end
    
  % norm_av_tool=av_rt_tool;

   opttool=fitoptions(eq, 'Lower', [0 0.25 max(norm_av_tool) min(norm_av_tool)], 'Upper', [Inf 2.25 max(norm_av_tool) min(norm_av_tool)], 'StartPoint', [Inf 1.25 max(norm_av_tool)/2  max(norm_av_tool)/2]);

   
try 
    
[fitobj_tool, goodness_tool, output_tool, convmsg_tool]=fit(distances,norm_av_tool',eq,opttool);
   rsq_tool=goodness_tool.rsquare;
   pps_border_tool=fitobj_tool.x0;
   
catch 

   rsq_tool=NaN;
   pps_border_tool=NaN;

end


% Bootper torso


    %new sample with only those trials 
    scramble_data_torso=ppsdata_def{3};

    
    %now we finally calculate the means 
    for g=1:5 %for the 5 distances 
    eachdist=find(scramble_data_torso(:,1)==distances(g)); %find trials for that distance
    mytrials=scramble_data_torso(eachdist,2);
    av_rt_torso(g)=mean(mytrials,'omitnan');
    std_rt_torso(g)=std(mytrials,'omitnan');
    zeta=[]; 
    zeta=abs((mytrials-av_rt_torso(g))/std_rt_torso(g));
    mytrials(:,3)=zeta>=cutoff; %writes 1 on the trials to discard 
    trialskeep=find(mytrials(:,3)==0);
    mytrials=mytrials(trialskeep,1);
    howmanyout_torso(g,:)=sum(zeta>cutoff);
    av_rt_torso(g)=mean(mytrials,'omitnan'); %new average for each distance of the cleaned data 
    std_rt_torso(g)=std(mytrials,'omitnan');
    end
    
    
    %we normalize the new averages between 0 and 1
   for l=1:length(av_rt_torso)
   norm_av_torso(l)=(av_rt_torso(l)-min(av_rt_torso))/(max(av_rt_torso)-min(av_rt_torso));
   end
    
   
     %norm_av_torso=av_rt_torso;
      opttorso=fitoptions(eq, 'Lower', [0 0.25 max(norm_av_torso) min(norm_av_torso)], 'Upper', [Inf 2.25 max(norm_av_torso) min(norm_av_torso)], 'StartPoint', [Inf 1.25 max(norm_av_torso)/2  max(norm_av_torso)/2]);

      

try 
    
   [fitobj_torso, goodness_torso, output_torso, convmsg_torso]=fit(distances,norm_av_torso',eq, opttorso);
   rsq_torso=goodness_torso.rsquare;
   pps_border_torso=fitobj_torso.x0;
   
catch 

   rsq_torso=NaN;
   pps_border_torso=NaN;

end
    


% Bootper whole


    %new sample with only those trials 
    scramble_data_whole=ppsdata_def{4};

    
    %now we finally calculate the means 
    for g=1:5 %for the 5 distances 
    eachdist=find(scramble_data_whole(:,1)==distances(g)); %find trials for that distance
    mytrials=scramble_data_whole(eachdist,2);
    av_rt_whole(g)=mean(mytrials,'omitnan');
    std_rt_whole(g)=std(mytrials,'omitnan');
    zeta=[]; 
    zeta=abs((mytrials-av_rt_whole(g))/std_rt_whole(g));
    mytrials(:,3)=zeta>=cutoff; %writes 1 on the trials to discard 
    trialskeep=find(mytrials(:,3)==0);
    mytrials=mytrials(trialskeep,1);
    howmanyout_whole(g,:)=sum(zeta>cutoff);
    av_rt_whole(g)=mean(mytrials,'omitnan'); %new average for each distance of the cleaned data 
    std_rt_whole(g)=std(mytrials,'omitnan');
    end
    
    
    %we normalize the new averages between 0 and 1
   for l=1:length(av_rt_whole)
   norm_av_whole(l)=(av_rt_whole(l)-min(av_rt_whole))/(max(av_rt_whole)-min(av_rt_whole));
   end
    

%  norm_av_whole=av_rt_whole;
   optwhole=fitoptions(eq, 'Lower', [0 0.25 max(norm_av_whole) min(norm_av_whole)], 'Upper', [Inf 2.25 max(norm_av_whole) min(norm_av_whole)], 'StartPoint', [Inf 1.25  max(norm_av_whole)/2  max(norm_av_whole)/2]);

   
try 
    
   [fitobj_whole, goodness_whole, output_whole, convmsg_whole]=fit(distances,norm_av_whole',eq, optwhole);
   rsq_whole=goodness_whole.rsquare;
   pps_border_whole=fitobj_whole.x0;
   
catch 

   rsq_whole=NaN;
   pps_border_whole=NaN;

end


figure(1)
hold on  

% --- SCATTER COLORS ---
scatter(distances, norm_av_pre,   80,  'MarkerEdgeColor',[1 0 0], 'MarkerEdgeAlpha',0.7,'LineWidth',2.5); % PRE = red
scatter(distances, norm_av_tool,  100, 'MarkerEdgeColor',[0 0 1], 'MarkerEdgeAlpha',0.7,'LineWidth',2.5); % PULL = blue
scatter(distances, norm_av_torso, 200, 'MarkerEdgeColor',[0 1 0], 'MarkerEdgeAlpha',0.7,'LineWidth',2.5); % TORSO = green
scatter(distances, norm_av_whole, 300, 'MarkerEdgeColor',[1 0 1], 'MarkerEdgeAlpha',0.7,'LineWidth',2.5); % WHOLE = magenta

% --- FITTED CURVES ---
pre   = plot(fitobj_pre,   'r'); pre.LineWidth   = 3;
tool  = plot(fitobj_tool,  'b'); tool.LineWidth  = 3;
torso = plot(fitobj_torso, 'g'); torso.LineWidth = 3;
whole = plot(fitobj_whole, 'm'); whole.LineWidth = 3;

xlim([0.22 2.28]);
ylim([0 1]);

xticks(0.25:0.5:2.25);
yticks(0:0.25:1);
yticklabels({'0.00','0.25','0.50','0.75','1.00'});

set(gca,'FontSize',20);
xlabel('Distance (m)');
ylabel('RTs (norm)');

set(gca,'TickLength',[0.03,0.01]);
set(gcf,'Color','w');
set(gca,'Color','w','xColor','k','YColor','k');
pbaspect([1 1 1]);
box off
set(gca,'TickDir','out');
set(gca,'LineWidth',3);

% --- ARROWS ---
arrowpre   = annotation('doublearrow','Color','r','HeadStyle','plain','LineWidth',1);
arrowpre.Parent = gca;
arrowpre.Position = [pps_border_pre, 0, 0, 0.15];

arrowtool  = annotation('doublearrow','Color','b','HeadStyle','plain','LineWidth',1);
arrowtool.Parent = gca;
arrowtool.Position = [pps_border_tool, 0, 0, 0.15];

arrowtorso = annotation('doublearrow','Color','g','HeadStyle','plain','LineWidth',1);
arrowtorso.Parent = gca;
arrowtorso.Position = [pps_border_torso, 0, 0, 0.15];

arrowwhole = annotation('doublearrow','Color','m','HeadStyle','plain','LineWidth',1);
arrowwhole.Parent = gca;
arrowwhole.Position = [pps_border_whole, 0, 0, 0.15];

legend('BASELINE','PULL','TORSO MOVEMENT','WHOLE BODY MOVEMENT','Location','best')

save('workspace_PPS.mat')

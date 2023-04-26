%% load ADHD and Healthy GRANDS
%[fileNames, pathName]=Z_getSetsFileNames;
clear all
cd 'D:\OneDrive\PhD Zangen\ADHD\NEW 2019 DRAFT'
load('D:\OneDrive\PhD Zangen\ADHD\NEW 2019 DRAFT\ADHD_stats_new2019.mat');% 'D:\OneDrive\PhD Zangen\ADHD\avi''s paper\BIG_FINAL_DATA_table.mat'
Final4=Final4(:,[1:23 26:27 30:34 46:40 42:end]);
%pathName='D:\WORKINGSET_D\ADHD\';
pathName='D:\OneDrive\PhD Zangen\ADHD\EEGLAB GRANDS BACKUP\TEP - 15ms\';
fileNames=[{'ADHD GRAND Merged datasets Unstich-basecorr.set'},...
     {'Healthy GRAND Merged datasets Unstich-basecorr.set'}];
ALLEEG=pop_loadset('filepath' , pathName, 'filename',fileNames);
% fileNames=[{'ADHD GRAND Merged datasets Unstich-basecorr.set'}];
% EEG = pop_loadset('filepath' , pathName, 'filename',fileNames);
% fileNames=[{'Healthy GRAND Merged datasets Unstich-basecorr.set'}];
% EEG2 = pop_loadset('filepath' , pathName, 'filename',fileNames);
% ALLEEG(1)=EEG; ALLEEG(2)=EEG2;
ALLEEG(1).condition = 'ADHD';
ALLEEG(2).condition = 'Healthy';
clear -vars pathName fileNames EEG EEG2

Final4.ADHD_smpts_total(Final4.group=='Healthy' & Final4.ADHD_smpts_total>70)=nan;
%% Sara-T iTBS pulse dosing Load GRANDS
clear all
cd 'D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\SP_preprocessed'
pathName='D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\SP_preprocessed';
%fileNames=[{'pre_theta_dose_GrandAverage.set'},...
%     {'post_theta_dose_GrandAverage.set'}];
fileNames=[{'pre600_theta_dose_GrandAverage.set'},...
     {'post600_theta_dose_GrandAverage.set'}];
ALLEEG=pop_loadset('filepath' , pathName, 'filename',fileNames);
%% load MST trial GRANDS

clear all
cd 'D:\WORKINGset_D\MST'
pathName='D:\WORKINGset_D\MST\';
fileNames=[{'pre_MST_MDD_SP_grand_GrandAverage.set'}, ...
    {'post_MST_MDD_SP_grand_GrandAverage.set'}];
ALLEEG = pop_loadset('filepath' , pathName, 'filename',fileNames);
ALLEEG(1).condition = 'PRE_MST';
ALLEEG(2).condition = 'POST_MST';
clear -vars pathName fileNames
%% Load Wellcome-LEAP GRANDS
%[fileNames, pathName]=Z_getSetsFileNames;
clear all
cd 'D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\ICA2\remcomp'
cd 'D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\ICA2\remcomp'

load('D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\ICA2\remcomp');% 'D:\OneDrive\PhD Zangen\ADHD\avi''s paper\BIG_FINAL_DATA_table.mat'
%Final4=Final4(:,[1:23 26:27 30:34 46:40 42:end]);
%pathName='D:\WORKINGSET_D\ADHD\';
pathName='D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\ICA2\remcomp';
fileNames=[{'WL_SP_PRE_GrandAverage.set'},...
     {'WL_SP_POST_GrandAverage.set'}];
ALLEEG=pop_loadset('filepath' , pathName, 'filename',fileNames);
ALLEEG(1).condition = 'PRE';
ALLEEG(2).condition = 'POST';
clear -vars pathName fileNames EEG EEG2


%% Group compare Historgrams 
co=brewermap(4,'PuOr'); %PuOr RdBu
ADHDcolor=co(1,:);Healthycolor=co(4,:); 
tabl=Final4;
vars=[16 17 18]; %16 17 18 26 37
for var=vars %  var=16
a=tabl{tabl.group=='Healthy',var}; a=a(isfinite(a)); b=tabl{tabl.group=='ADHD',var}; b=b(isfinite(b)); 
% a = bootstrp(1000,@(x)[mean(x) std(x)],a)
%a = bootstrp(100,@nanmean,a); b=bootstrp(100,@nanmean,b);
%pd = makedist('Normal','mu',nanmean(a),'sigma',nanstd(a))
[fia,xia] = ksdensity(a);[fib,xib] = ksdensity(b);
figure; hold on;
%plot(xia,fia); plot(xib,fib);
%hist(a) gca.Healthycolor; hist(b,ADHDcolor)
fill(xia,fia,Healthycolor,'EdgeColor',[Healthycolor.*0.6],'FaceVertexAlphaData',0.1,'FaceAlpha',0.55); 
fill(xib,fib,ADHDcolor,'EdgeColor',[ADHDcolor.*0.6],'FaceVertexAlphaData',0.1,'FaceAlpha',0.55);
legend('Healthy','ADHD'); legend boxoff; title(strrep(tabl.Properties.VariableNames(var),'_',' '));
hold off
%a = datasample(a,1000); b = datasample(b,1000)
%binsa=ceil(sqrt(sum(isfinite(a)))); binsb=ceil(sqrt(sum(isfinite(b))));
%figure; hold on;
%generatePDF(a,Healthycolor,'area',binsa); generatePDF(b,ADHDcolor,'area',binsb)
%legend('Healthy','ADHD'); legend boxoff; title(strrep(tabl.Properties.VariableNames(var),'_',' '));
%hold off;
end
%% Group descriptive statistics
var=26
tabl=Final5;
[h, p, ci, stats_tt]=ttest2(tabl{tabl.group=='Healthy',var},tabl{tabl.group=='ADHD',var})
sta = varfun(@(x) [mean(x,'omitnan'); std(x,'omitnan')],tabl,'GroupingVariables','group' ,'InputVariables',var );
sta.Properties.RowNames={'Mean', 'STD','Mean2', 'STD2'}
%
var=37;
tabl=Final5;
sum(isfinite(tabl{tabl.group=='Healthy',var}))
sum(isfinite(tabl{tabl.group=='ADHD',var}))
%% POP_COMPARE EEGLAB

time= [-50 400]
pop_comperp( ALLEEG, 1, [1 2] ,[],'addavg','off','addstd','off','addall','on','diffavg','off','diffstd','off','tplotopt',{'ydir' 1},'tlim',time );
%  print(gcf,'-dpdf','-painters','-bestfit','-r600','channel_array_pre_post_compare-50_100');
%% subject remove
subrem=256
posit=find(strcmp(ALLEEG(1).subjects,string(subrem)))
ALLEEG(1).data(:,:,posit) = [];
ALLEEG(1).epoch(:,posit) = [];
ALLEEG(1).subjects(:,posit) = [];
Final4=Final4(~(Final4.Subjects==subrem),:);

%
isoutlier(Final4.ADHD_smpts_total(Final4.group=='Healthy'))
Final4.ADHD_smpts_total(Final4.group=='Healthy' & Final4.ADHD_smpts_total>70)=nan;
isoutlier(Final4.ADHD_index(Final4.group=='Healthy'))
Final4.ADHD_index(Final4.group=='Healthy' & Final4.ADHD_index>70)=nan;
isoutlier(Final4.TOTAL(Final4.group=='Healthy'))
Final4.TOTAL(Final4.group=='Healthy' & Final4.TOTAL>70)=nan;
%% subject remove
subrem=256
posit=find(strcmp(ALLEEG(1).subjects,string(subrem)))
EEG.data(:,:,posit) = [];
EEG.epoch(:,posit) = [];
EEG.subjects(:,posit) = [];
%% Remove channel from subject
% figure; topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', ALLEEG(1).chaninfo);
      
temp=EEG
temp.data=[]
temp.data(:,:,1)=EEG.data(:,:,41)
temp.data(:,:,2)=EEG.data(:,:,41)
temp=pop_select(temp, 'trial', [1 2]);
temp=pop_select(temp,'nochannel', 1);
plot(temp.data(:,:,1)')

elecnum=elecName(EEG,{'po8' 'po4'});
EEG=pop_select(EEG,'nochannel', elecnum);
load('D:\Itay-Google-Drive\MATLAB\LAB_MatlabScripts\chanlocations.mat')
EEG = pop_interp(EEG, ChanLocations, 'spherical');

    %export_fig Headplot_electrodes_N3.pdf -q101 -painters
    
ALLEEG(1) = pop_select( ALLEEG(1),'nochannel',[13 19]);
ALLEEG(2) = pop_select( ALLEEG(2),'nochannel',[13 19])
%% trim and interpulate data 

trim=[-2 20]; %DoublePulseINT=0
ALLEEG(1)= Triminterp(ALLEEG(1),trim);
%% Ratio between 2 table variable
Final.TEP30_TEP180_ratio=table2array(Final(:,27))./table2array(Final(:,28))
%% MST trial
clearvars -except Final2 Final4 fileNames pathName ALLCOM ALLEEG CURRENTSET CURRENTSTUDY chanlocs lo hi EEG LASTCOM PLUGINLIST STUDY eeglabUpdater CARERP CAARS Final Final2
dataset=[1 2]
%ALLEEG(1).condition = 'ADHD';%      ALLEEG(1).condition = 'POST';
%ALLEEG(2).condition = 'Healthy'; %     ALLEEG(2).condition = 'PRE';
Groups= {ALLEEG(dataset(1)).condition ALLEEG(dataset(2)).condition}%{'Pre' 'Post'};
wavgrph='on'
violinplt='on'
fixfit='o'
facelift='o'
LMFPgrph='on'
GMFPgrph='o'
ISP='o'
bargrph='o'
rectyfied='o'
stichsize=17
% P30: 20-40ms, N45: 40-50ms, P60: 50-70ms, N100: 100-130ms, P180: 160-250ms 
name='P30'
statwin=[20 90]; % for N200 [160 220]

% graamplot.export('file_name',strrep(['violin ' name ' ' num2str(statwin)],' ','_'),'file_type','pdf')
elec={'fc3' 'fc4' 'f3' 'f5'} % 'f2' 'fc2' 'fc4' 'fc6' 'fcz' 'fc2' 'fc4' 'fz' 'f2' 'f4' 'fc2' 'fc4' 'fcz'  'fcz' 'fz'     for N200  {'fc2','fc4' ,'f2' , 'f4'} {'f2' 'f4' 'fz'  'f4' 'fc4'      'fcz' 'fc2' 'fc4'      'f2' 'fz' 'f2' 'f4' 'fz'  'fcz' 'fc2' 'fc4'
timeWin=[-100 500];
amp=20
method='LMFP'; %'mean''admean''AUC''LMFP'
rem_outlie='off'
peaktype='positive' %positive
smallWinHalfLength=2

[stats, head, TEPVec1, graamplot]=TEPstats(ALLEEG, dataset, Groups,name, fixfit, facelift, rem_outlie, violinplt, ISP,rectyfied, wavgrph,LMFPgrph,GMFPgrph, bargrph, timeWin, statwin,stichsize, elec, amp, method, peaktype, smallWinHalfLength);
%% LMFP ADHD cross sectional
clearvars -except Final2 Final4 fileNames pathName ALLCOM ALLEEG CURRENTSET CURRENTSTUDY chanlocs lo hi EEG LASTCOM PLUGINLIST STUDY eeglabUpdater CARERP CAARS Final Final2
dataset=[1 2]; Groups={ALLEEG(dataset(1)).condition ALLEEG(dataset(2)).condition};%{'Pre' 'Post'};
wavgrph='o'; fixfit='o'; facelift='off';stichsize=17;
LMFPgrph='on'; GMFPgrph='o'; ISP='o'; rectyfied='o';
violinplt='on'; bargrph='o';
name='P30'; %P30: 20-40ms, N45: 40-50ms, P60: 50-70ms, N100: 100-130ms, P180: 160-250ms 
statwin=[25 45]; % for N200 [160 220]
elec={'fc6' 'f6' 'ft8' 'f8' } ;%   'ft8' 'f8' 'fc4' 'fc6' 'f4' 'f6' 'fc2' 'fc4' 'fc6' 'f4' 'f6' 'f2' 'fc2' 'fc4' 'fc6' 'fcz' 'fc2' 'fc4' 'fz' 'f2' 'f4' 'fc2' 'fc4' 'fcz'  'fcz' 'fz'     for N200  {'fc2','fc4' ,'f2' , 'f4'} {'f2' 'f4' 'fz'  'f4' 'fc4'      'fcz' 'fc2' 'fc4'      'f2' 'fz' 'f2' 'f4' 'fz'  'fcz' 'fc2' 'fc4'
% figure; topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', ALLEEG(1).chaninfo);
timeWin=[-100 400]; amp=10;
method='LMFP'; %'mean / std / rectified / AUC / GMFP / LMFP / admean
peaktype='positive'; smallWinHalfLength=2; %positive
rem_outlie='on';
[stats, head, TEPVec1, graamplot]=TEPstats(ALLEEG, dataset, Groups,name, fixfit, facelift, rem_outlie, violinplt, ISP,rectyfied, wavgrph,LMFPgrph,GMFPgrph, bargrph, timeWin, statwin,stichsize, elec, amp, method, peaktype, smallWinHalfLength);

%chanlocs=load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\chanlocs60');
% graamplot.export('file_name',strrep(['violin LMFP' name ' ' num2str(statwin)],' ','_'),'file_type','pdf')
%load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat'); chanlocs=chanlocs62; %for the sara T THETA-burst analysis
%lo=-4;hi=4;[stat] = clusterSTATS(head{1}, head{2},lo,hi,chanlocs,stats)
%Final2=stats;
if ~exist('Final2')
    Final2=stats;
    Final5=outerjoin(Final4,Final2,'MergeKeys',true,'Keys',{'Subjects' 'group'});
else
    Final2=outerjoin(Final2,stats,'MergeKeys',true,'Keys',{'Subjects' 'group'});
    Final5=outerjoin(Final4,Final2,'MergeKeys',true,'Keys',{'Subjects' 'group'});
end
 %  figfig(strrep(['ADHD_Control_LMFP_timecourse ' num2str(statwin)],' ','_'))
%  BIG_Final=outerjoin(BIG_Final,Final2(:,[1 3]),'MergeKeys',true);
%Final2{:,13}=str2num(str2mat(string(Final2{:,1})))
% Final3=outerjoin(Final,Final2(:,[1 2 8:12]),'MergeKeys',true);
%ttest2(Final2(Final2.Group==ADHD),Final2(Final2.Group==Healthy)
%% Sara-Manon iTBS pulse dosing - TEP cpmponents 
clearvars -except Final2 Final4 fileNames pathName ALLCOM ALLEEG CURRENTSET CURRENTSTUDY chanlocs lo hi EEG LASTCOM PLUGINLIST STUDY eeglabUpdater CARERP CAARS Final Final2
dataset=[1 2]
Groups= {ALLEEG(dataset(1)).condition ALLEEG(dataset(2)).condition}%{'Pre' 'Post'};
wavgrph='on'; violinplt='on'; fixfit='o'; facelift='o';LMFPgrph='o'; GMFPgrph='o'; ISP='o'
bargrph='o';rectyfied='o'; stichsize=17
name='P200'
statwin=[170 220]; % for N200 [160 220] P30 (20 35), N45 (40 50), P60 (55 65), N100 (85 125) and P200 (170 220). 
% graamplot.export('file_name',strrep(['violin ' name ' ' num2str(statwin)],' ','_'),'file_type','pdf')
elec=   { 'AF3' 'F1' 'F3' 'FC1' 'FC3'} %  elec= {'AF3' 'F1' 'F3' 'FC1' 'f4' 'FC3' 'fc2','fc4' ,'f2' , 'f4' 'f4'}
timeWin=[-100 300]; amp=8
method='mean'; %'mean''admean''AUC''LMFP'
rem_outlie='off'; peaktype='positive'; smallWinHalfLength=2

[stats, head, TEPVec1, graamplot]=TEPstats(ALLEEG, dataset, Groups,name, fixfit, facelift, rem_outlie, violinplt, ISP,rectyfied, wavgrph,LMFPgrph,GMFPgrph, bargrph, timeWin, statwin,stichsize, elec, amp, method, peaktype, smallWinHalfLength);
% figfig('pre-post_600pulse')
%chanlocs=load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\chanlocs60');

% graamplot.export('file_name',strrep(['violin ' name ' ' num2str(statwin)],' ','_'),'file_type','pdf')

%load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat'); chanlocs=chanlocs62; %for the sara T THETA-burst analysis
%lo=-4;hi=4;[stat] = clusterSTATS(head{1}, head{2},lo,hi,chanlocs,stats)
%Final2=stats;

%Final2{:,13}=str2num(str2mat(string(Final2{:,1})))
% Final3=outerjoin(Final,Final2(:,[1 2 8:12]),'MergeKeys',true);
%ttest2(Final2(Final2.Group==ADHD),Final2(Final2.Group==Healthy)
%export_fig TEP_ECT_PRE_POST.pdf -q101 -painters -append -transparent
%figfig('GMFP_Theta-Dose_PRE_POST.pdf')
%% ttest between groups in spesific variable
var=size(Final2,2)
var=[8:12]
Final2(41,var)={nan};
Final2(5,var)={nan};
% Final2=Final4 % var=26

[h p ci stats]=ttest2(table2array(Final2(Final2.group=='ADHD',var)),table2array(Final2(Final2.group=='Healthy',var)))
[h p_cor]=fdr(p,0.05)
[h p_cor p_adj]=TMS_fdr(p,0.05)
hist(table2array(Final2(Final2.group=='Healthy',16)))

var=4%size(ERSPstats,2)
[h p]=ttest2(table2array(ERSPstats(ERSPstats.group=='ADHD',var)),table2array(ERSPstats(ERSPstats.group=='Healthy',var)))
hist(table2array(ERSPstats(ERSPstats.group=='Healthy',var)))
hist(table2array(ERSPstats(ERSPstats.group=='ADHD',var)))
ERSPstats(30,2:5)={nan}
% export_fig ADHD_HEALTHY_TEP_P180_bar_basecorr.PNG -q101
% print(gcf, '-dpdf','-painters','-bestfit',  'TEP_ECT_PRE_POST.pdf')
%figure; hold on; scatter(repmat(1,size(Final2(Final2.group=='ADHD',var),1),1),table2array(Final2(Final2.group=='ADHD',var)),'r')
%scatter(repmat(2,size(Final2(Final2.group=='Healthy',var),1),1),table2array(Final2(Final2.group=='Healthy',var)),'b')
%  writetable(Final2,'ADHD_cont_TEP.xlsx')
mean(Final2{Final2.group=='ADHD',8})
std(Final2{Final2.group=='ADHD',8})
mean(Final2{Final2.group=='Healthy',8})
std(Final2{Final2.group=='Healthy',8})
%% Linear regression Model 
vars=[6 24]; %16
resp=10;
tabl=table; tabl=Final5;
filtnan=[]; predictors=[]; response=[]; model=[]; mdl=[];
filtnan=all(isfinite(tabl{:,vars}),2);
predictors=tabl(filtnan,[vars]);
if strcmpi(class(tabl{filtnan,resp}),'categorical') || ischar(tabl{filtnan,resp})
 response=table(grp2idx(tabl{filtnan,resp}),'VariableNames',tabl.Properties.VariableNames(resp));
else %isfinite(Final4{filtnan,resp})
 response=tabl(filtnan,resp);   
end
model=char(join([tabl.Properties.VariableNames{resp} '~' join({tabl.Properties.VariableNames{vars}},'+')]));
mdl = fitlm([predictors, response],model)
%mdl = fitglm([predictors, response],model)
 
stats=anova(mdl,'summary')

% one way anova
[p,tbl,stats]=anova1(predictors{:,:}, tabl{filtnan,resp},'summary')
    [c,m,h,nms] = multcompare(stats,'Estimate','column')
    
    
 %   [p,tbl,stats]=anova2(predictors{:,:}, 91,'summary')
%% repeated ANOVA - not necessery here
vars=[16 50]
resp=2;
filtnan=all(isfinite(Final4{:,vars}),2);
predictors=Final4(filtnan,[vars]);
response=Final4(filtnan,resp);
model=[Final4.Properties.VariableNames{vars(1)} '-' Final4.Properties.VariableNames{vars(2)} '~' Final4.Properties.VariableNames{2}];
% repeated=vars;
rm = fitrm([predictors response], model, 'WithinDesign', [1 2]);
no = ranova(rm)
post = multcompare(rm,'group')

% another option 
aa=readtable('A:\WorkingSet\ThetaBurst_DOSE_ST\test.xlsx')
aaa=aa(:,2:7);
WithinStructure = table(categorical([1 2 1 2 1 2]'),categorical([1 1 2 2 3 3]'),'VariableNames',{'PrePost','Dose'});
rm = fitrm(aaa, 'pre6,post6,pre12,post12,pre18,post18~1','WithinDesign',WithinStructure);
ranovatable = ranova(rm,'WithinModel','PrePost*Dose');
plotprofile(rm,'Dose')
post = multcompare(rm,'Dose')
boxplot(aaa{:,:})
%% removing subjects with NANs
Final3(sum(isnan(table2array(Final3(:,2:end))),2)>=1,:)=[]
%% correlation matrix
correl=outerjoin(Final(:,[1:2,10:18]),stats,'MergeKeys',true);
%correl=Final
%[c p]=corr(table2array(Final(:,vars(1))),table2array(Final(:,vars(2))))
%[cormatADHD.R cormatADHD.P)=corrcoef(table2array(corrsADHD),'rows','pairwise') %{Final.Properties.VariableNames(relevant) num2cell(cormat)}
%[cormatHealthy.R cormatHealthy.P)=corrcoef(table2array(corrsHealthy),'rows','pairwise')
[cormat.R cormat.P]=corrcoef(table2array(correl(:,3:end)),'rows','pairwise');
corMat=nan(size(cormat.R));
corMat(cormat.P<0.05)=cormat.R(cormat.P<0.05);
corMat=array2table(corMat,'VariableNames',correl.Properties.VariableNames(3:end),'RowNames',correl.Properties.VariableNames(3:end)')
%% 2-groups correlations
vars=[26 48];
alpha=0.05;
disregard_zeros='on';
disregard_num=0; %should be 0 to disable
rem_outlie='on';
cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
%[ 
[corMatADHDF, corMatHealthyF, corMatF] = Corr_stats(Final2,vars, 'Healthy', 'MDD', alpha, cortype,disregard_zeros, disregard_num, rem_outlie );
% Final3{:,30:34}=double(Final3{:,30:34})    Final3{:,34}
%% Correlations -2-groups ADHD Controls 
%Final4=Final2
if ~exist('Final2')
    Final5=Final4;
   
elseif exist('Final2') && exist('Final4')
    Final5=outerjoin(Final4,Final2,'MergeKeys',true,'Keys',{'Subjects' 'group'});
end
% 
varcorrect=0; var=10;
if varcorrect==0

%Final4{:,varco+1}=10*log10(Final4{:,varco}); %LOG a VAR
lm=fitlm(Final5.Gender,Final5{:,var}); %partialling gender
Final5{:,var}=lm.Residuals.Raw;
lm=fitlm(Final5.Age,Final5{:,var}); %partialling gender
Final5{:,var}=lm.Residuals.Raw;
end
%Final4{Final4{:,16}>=85,16}=nan
%Final5{:,end+1}=abs(10*log10(Final5{:,38}))
%Final4{:,50}=Final4{:,16}; Final4{Final4.group=='ADHD' & Final4{:,16}<60 ,50}=nan 
vars=[6 17]; %16 39    6 24   10 24    10 6  16 37    5 26
alpha=0.1;
disregard_zeros='o';
disregard_num=0; %should be 0 to disable
rem_outlie='on';
cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
%Final4=outerjoin(Final2,Final,'MergeKeys',true);
[corMatHealthy, corMatADHD, corMat, Figu] = Corr_stats(Final5,vars, 'Healthy', 'ADHD', alpha, cortype,disregard_zeros, disregard_num, rem_outlie );
figfig('supp_corr')
% compare_correlation_coefficients(-0.53,-0.56,17,11)
% final3(sum(~isnan(table2array(Final3(:,2:end))),2)>1,:)
% badvar=67;     thresh=7
% Final2(table2array(Final2(:,badvar))<thresh,badvar)={nan};       %OUTLIERS
% Final(table2array(Final(:,badvar))>thresh,badvar)={nan};       %OUTLIERS
% writetable(correlations_table,'correlations_table.xlsx')
%% Correlation-  one group 
vars=[15 16];
alpha=0.05;
cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
disregard_zeros='on';
disregard_num=90; %should be 0 to disable
rem_outlie='on';

[correlations]=Corr_simp(Final4, vars,alpha,cortype,disregard_zeros,disregard_num, rem_outlie);
% export_fig MDD_rTMS_baseline_correlations_sLORETA-GMFP.PDF -append -q101 -painters
% print(gcf, '-dpdf','-painters','-bestfit',  'TEP_ECT_PRE_POST.pdf')
%% Non Linear regression model

if ~exist('Final2')
    Final4=Final2;
   % Final2=stats;
elseif exist('Final2') & exist('Final4')
    Final4=outerjoin(Final4,Final2,'MergeKeys',true,'Keys',{'Subjects' 'group'});
end

outlie='o';
covar='o';
for aaa=[16]
for ppp=[37]
VARS=[aaa ppp];
%Final4{Final4{:,VARS(1)}>=200,VARS(1)}=nan;
%Final4{Final4{:,VARS(2)}>=6,VARS(2)}=nan;

corrs=Final4(:,VARS);

%corrs{:,2}=10.*log10(corrs{:,2});
%corrs{:,1}=10*log10(corrs{:,1});


if strcmp(covar,'on')
lm=fitlm(Final4.Gender,corrs{:,1}); %partialling gender
corrs{:,1}=lm.Residuals.Raw;

lm=fitlm(Final4.Gender,corrs{:,2}); %partialling gender
corrs{:,2}=lm.Residuals.Raw;

lm=fitlm(Final4.Age,corrs{:,1}); %partialling gender
corrs{:,1}=lm.Residuals.Raw;

lm=fitlm(Final4.Age,corrs{:,2}); %partialling gender
corrs{:,2}=lm.Residuals.Raw;
end

if strcmp(outlie,'on')
corrs{isoutlier(corrs{:,1}),1}=nan;
corrs{isoutlier(corrs{:,2}),2}=nan;
end

[X, Y] = prepareCurveData( corrs{:,1},corrs{:,2});

%modelfun=@(a,x)a(1)*x.^a(2) % power function
%modelfun=@(a,x)a(1)*x+a(2) % linear function - didn't try
modelfun=@(a,x) a(1).*x.^2+a(2).*x+a(3); % polynomal 2nd order function
mdl=fitnlm(X,Y,modelfun,[0; 0; 0]);
%modelout=mdl.Residuals.Standardized<-3 | mdl.Residuals.Standardized>3;
%corrs{modelout,:}=NaN;
%[p F]= coefTest(mdl);

ssr = max(mdl.SST - mdl.SSE,0);
dfr = mdl.NumEstimatedCoefficients - 1;
dfe = mdl.NumObservations - 1 - dfr;
f = (ssr./dfr) / (mdl.SSE/dfe);
p = fcdf(1./f,dfe,dfr); % upper tail
if p<0.1
%h=plotSlice(mdl);
% Xaxis=[nanmin(X(1,:)):0.1:nanmax(X(1,:))]';
%Xaxis=[nanmin(corrs{:,1}):0.1:nanmax(corrs{:,1})]';
Xaxis=[nanmin(X(:,1)):0.1:nanmax(X(:,1))]';
[pfun, pCI]=predict(mdl, Xaxis); 
figure;  hold on,
co=brewermap(4,'PuOr'); %PuOr RdBu
ADHDcolor=co(2,:);
ADHDcolor1=co(1,:);%'r'; %[0 0 0];
Healthycolor=co(3,:); % 'b' ;%[166/255 166/255 166/255];
Healthycolor1=co(4,:);
a=scatter(corrs{Final4.group=='ADHD',1},corrs{Final4.group=='ADHD',2},85,'o','LineWidth',0.75,'MarkerFaceColor',ADHDcolor,'Markerfacealpha',.85,'MarkerEdgeColor',ADHDcolor1);
h=scatter(corrs{Final4.group=='Healthy',1},corrs{Final4.group=='Healthy',2},85,'o','LineWidth',0.75,'MarkerFaceColor',Healthycolor,'Markerfacealpha',.6,'MarkerEdgeColor',Healthycolor1);
h2=plot(Xaxis,pfun,'--','color', [0.1 0.1 0.1 0.8], 'LineWidth', 2);
plot(Xaxis,pCI,'--','Color',[0.1 0.1 0.1 0.5] )
xlabel(strrep(Final4.Properties.VariableNames{VARS(1)},'_',' '));
ylabel(strrep(Final4.Properties.VariableNames{VARS(2)},'_',' '));
legend([a h h2],['ADHD'],['Healthy'],['ploynomial trend' ' (n=' num2str(mdl.NumObservations) '), R='... 
            num2str(mdl.Rsquared.Ordinary.^0.5, 2) ', p=' num2str(p,2)] );
hold off
%export_fig PolyREG_ADHD_CONT.pdf -q101 -painters -append %
%close all
end
end
end

% figure; 
% plot(mdl)
% figure;
% plotDiagnostics(mdl,'contour')
% plotResiduals(mdl)
%% ROC curve 
%Final4=BIG_Final(:,[1:23 26:27 30:34 46:40 42:end]);

% discr=16;
% Final5=table;
% Final5=Final4(Final4.group=='ADHD',:);
% Final5{:,49}=nan(size(Final5,1),1);
% Final5{Final5{:,discr}> prctile(Final5{:,discr},75),49}=true;
% Final5{Final5{:,discr}< prctile(Final5{:,discr},75),49}=false;
VAR_corr=0;
resp=[2];
vars=[6 37]; % [6 7 8 37]
%[6 7 8 24 41] 4 6 16 21 26 6 21 48 16 19 26 29 30 %  18 6 26 48

Final6=table; Final6=Final5;
Final6=Final6(~sum(ismissing(Final6{:,[resp]}),2)>0,:);
Final6=Final6(~sum(ismissing(Final6{:,[vars]}),2)>0,:);
if VAR_corr==1
for jojo=1:size(vars,2)
lm=fitlm(Final6.Gender,Final6{:,vars(jojo)}); %partialling gender
Final6{:,vars(jojo)}=lm.Residuals.Raw;
lm=fitlm(Final6.Age,Final6{:,vars(jojo)}); %partialling gender
Final6{:,vars(jojo)}=lm.Residuals.Raw;
%lm=fitlm(Final5.Beck_Total,Final5{:,vars(jojo)}); %partialling gender
%Final5{:,vars(jojo)}=lm.Residuals.Raw;
end
end
% try out the fitcsvm function!!!
clsfun='linear'; % logistic linear
ROCcurve(Final6,vars,resp,clsfun)
%% Partial Correlation
VAR1=[6]; 
VAR2=[10];
contVAR=[24];
tabl=Final5;
natu=isfinite((sum(tabl{:,[VAR1 VAR2 contVAR]},2)));

%subs=natu & grp;

[rho pval]=partialcorr(tabl{[natu & grp],VAR1} ,tabl{[natu & grp],VAR2}, tabl{[natu & grp],contVAR})
%% 3D - Scatter plot 
tabl=Final4;
VAR1=[6]; 
VAR2=[40];
VAR3=[24];
Groups=unique(tabl.group);
grp1=[all(isfinite(tabl{:,[VAR1 VAR2 VAR3]}),2) & tabl.group==Groups(1)];
grp2=[all(isfinite(tabl{:,[VAR1 VAR2 VAR3]}),2) & tabl.group==Groups(2)];
co=brewermap(4,'PuOr'); %PuOr RdBu
ADHDcolor=co(1,:);  Healthycolor=co(4,:); 
ff=figure; hold on
a=scatter3(tabl{grp1,VAR1},tabl{grp1,VAR2},tabl{grp1,VAR3},...
    160,'o','LineWidth',0.1,'MarkerFaceColor',ADHDcolor,'Markerfacealpha',.4,'MarkerEdgeColor',[ADHDcolor.*0.7]);
b=scatter3(tabl{grp2,VAR1},tabl{grp2,VAR2},tabl{grp2,VAR3},...
    160,'o','LineWidth',0.1,'MarkerFaceColor',Healthycolor,'Markerfacealpha',.4,'MarkerEdgeColor',[Healthycolor.*0.7]);
xlabel(strrep(['Var ' num2str(VAR1) '- ' tabl.Properties.VariableNames{VAR1}],'_',' '));
ylabel(strrep(['Var ' num2str(VAR2) '- ' tabl.Properties.VariableNames{VAR2}],'_',' '));
zlabel(strrep(['Var ' num2str(VAR3) '- ' tabl.Properties.VariableNames{VAR3}],'_',' '));
set(gca,'CameraPosition',[-50 -440 -600],'FontName','Helvetica Neue','FontSize',10)
legend([a b], [char(Groups(1)) ' n=' num2str(sum(grp1))],[char(Groups(2)) ' n=' num2str(sum(grp2))],'FontSize',10);
ff.Renderer='painters';
hold off;
%% Feature selection - check this OUT!!!
vars=[16 6 26 48];  %4 7 8 26 30 48
resp=[2];
mdl=fscnca(Final4{:,vars},Final4{:,resp})%,...
%'FitMethod','exact','Solver','sgd'),...
 %   'Standardize',true);
figure; 
plot(mdl.FeatureWeights,'o','LineWidth',0.75,'MarkerSize',7,'MarkerFaceColor',[1,0.2,0.2],'MarkerEdgeColor',[1,0,0])

grid on
set(gca,'xtick',[1:size([3:size(Final4,2)],2)],'xticklabel',strrep(Final4.Properties.VariableNames(vars),'_',' '))
xtickangle(gca,45)
ylabel('Feature weight')

mdl=fsrftest(Final4{:,vars},Final4{:,resp})%,...
idx = fscmrmr(Tbl,ResponseVarName)
%% logistic classifier
var=[16 18 35] % Final2 table variables to classify
[B,dev,stats] = mnrfit(table2array(Final4(:,var)),Final2.group);
%% ratio between coloumn
Final2.TEP30_TEP180_ratio=table2array(Final2(:,27))./table2array(Final2(:,28))
%% creating SPL
headplot('setup',EEG.chanlocs,['GW64.spl'], 'chaninfo', EEG.chaninfo, 'meshfile','mheadnew.mat','transform',[-0.35579 -6.3369 12.3705 0.053324 0.018746 -1.5526 1.0637 0.98772 0.93269]);
% headplot('setup',ALLEEG(1).chanlocs,['GW64.spl'], 'chaninfo', ALLEEG(1).chaninfo, 'meshfile','mheadnew.mat','transform',[-0.35579 -6.3369 12.3705 0.053324 0.018746 -1.5526 1.0637 0.98772 0.93269]);
%% topoplot time series
co=brewermap(64,'OrRd');% http://colorbrewer2.org/       %   YlOrBr RdYlBu PuOr RdBu   OrRd   RdYlBu   RdBu
pop_topoplot(ALLEEG(1),1, [20:10:50] ,['timeseries-' ALLEEG(1).condition],[5 10] ,0,'electrodes','off','colormap',co, 'style', 'map', 'headrad', 0,'shading','interp','maplimits',[-3.5 3.5]);%,'maplimits',[0 0.1]
colormap(co)
pop_topoplot(ALLEEG(2),1, [20:10:50] ,['timeseries-' ALLEEG(2).condition],[5 10] ,0,'electrodes','off', 'style', 'map', 'headrad', 0,'shading','interp','colormap',co,'maplimits',[-3.5 3.5]);%,'maplimits',[0 0.1]
colormap(co)
% topoplot time series - 2 dataset subtration
EEGsub=ALLEEG(1); EEGsub.data=[]; EEGsub.data(:,:,1)=mean(ALLEEG(2).data,3)-mean(ALLEEG(1).data,3); EEGsub.data(:,:,2)=EEGsub.data(:,:,1); EEGsub.trials=2;
pop_topoplot(EEGsub,1, [20:10:50] ,['timeseries-' ALLEEG(2).condition],[5 10] ,0,'electrodes','off', 'style', 'map', 'headrad', 0,'shading','interp','colormap',co,'maplimits',[-3.5 3.5]);%,'maplimits',[0 0.1]
colormap(co)
%% stats headplot Ttests - Significance map - the "head" variable is produced in TEPstats function
    % clearvars -except fileNames pathName ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG LASTCOM PLUGINLIST STUDY eeglabUpdater CARERP CAARS Final Final2
       [~, p]=ttest2(head{1},head{2}, 'Dim', 2);
    p_log=-log10(p);
    figure; headplot(p, 'GW64.spl','maplimits',[0 0.1])
    figure; headplot(p_log, 'GW64.spl','maplimits',[-2 2],'electrodes','off');
    figure; topoplot(p_log, EEG.chanlocs,'maplimits',[-2 2]);colorbar
    % export_fig ADHD_HEALTHY_P180_sigHead.PNG -q101
    % export_fig ADHD_HEALTHY_P180_sigtopo_for_scale.PDF -q101 -painters
    % export_fig topoplot_ADHD_HEALTHY_120-300ms.PNG -q101 -painters
    % print(gcf, '-dpdf','-painters','-bestfit',  'TEP_ECT_PRE_POST.pdf')
    
    [~, p]=ttest2(squeeze(headminmax{1}),squeeze(headminmax{2}), 'Dim', 2);
    p_log=-log10(p);
    figure; headplot(p, 'GW64.spl','maplimits',[0 1])
    figure; headplot(p_log, 'GW64.spl','maplimits',[-2 2],'electrodes','off');
    figure; topoplot(p, EEG.chanlocs,'maplimits',[0 1],'style','map','shading','interp');colorbar ; title('P Value')
    figure; topoplot(p_log, EEG.chanlocs,'maplimits',[-2 2],'style','map','shading','interp');colorbar ; title('log10(P Value)')
    % export_fig ADHD_HEALTHY_MIN_MAX_log10_P_value.PNG -q101
    % print(gcf, '-dpdf','-painters','-bestfit',  'TEP_ECT_PRE_POST.pdf')
%% STD - standard deviation timeseries topoplot
clear EEGi
for gr=1:2
EEG=ALLEEG(gr);
%func= @(x) std(x,0,1);
data=zeros(size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));
LMFP=zeros(size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));
for j=1:size(EEG.data,3)
for i=1:size(EEG.data,1)
%data(i,:,j) =blockproc(EEG.data1(i,:,j),[1 6],func);
 data(i,:,j) = movstd(EEG.data(i,:,j),6);
 LMFP(i,:,j)= smoothdata(double(sqrt([EEG.data(i,:,j)-mean(EEG.data(i,:,j),2,'omitnan')]).^2),'gaussian',15);;
 %LMFP(i,:,j)=double(squeeze(nanmean(sqrt(nansum([t(elecnum,timeInd,:)-nanmean(t(elecnum,timeInd,:),2)].^2,1)./size(elecnum,2)),3)));
 %LMFP(:,:,i)=smoothdata(LMFP(:,:,i),'gaussian',15);
end
end
EEGi(gr)=EEG;
EEGi(gr).data=LMFP;
%ALLEEG(gr).data=data;
limits=[0 0.5];
%pop_topoplot(EEGi(gr),1, [16:10:120] ,['STD timeseries-' ALLEEG(gr).condition],[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0 );%,'maplimits',[0 0.1]
end

% LMFP-like topoplot
EEGi(3)=EEG; EEGi(3).data=[]; EEGi(3).trials=2;
size(EEGi(1).data)
EEGi(3).data=mean(EEGi(2).data,3)-mean(EEGi(1).data,3); EEGi(3).data(:,:,2)=mean(EEGi(2).data,3)-mean(EEGi(1).data,3);
gr=2
pop_topoplot(EEGi(gr),1, [20:10:60] ,['STD timeseries-' EEGi(gr).condition],[5 10] ,0,'electrodes','off', 'style', 'map', 'headrad', 0,'shading','interp','maplimits',[-3.5 3.5]);%,'maplimits',[0 0.1]
% figfig('timeseries_std_topoplots')
%% CSD (current source density) - trying to deduce the current from the given potentials(??)
% del2map
HeadCSD1=del2map(double(head{1}),EEG.chanlocs);
HeadCSD2=del2map(double(head{2}),EEG.chanlocs);

%HeadCSD=del2map(double(EEG.data),EEG.chanlocs)
figure; hold on; topoplot(mean(HeadCSD1,2), EEG.chanlocs);  title([ALLEEG(1).setname], 'interpreter', 'none'); colorbar; hold off;
figure; hold on; topoplot(mean(HeadCSD2,2), EEG.chanlocs); title([ALLEEG(2).setname], 'interpreter', 'none');colorbar; hold off;

figure; headplot(HeadCSD2, 'GW64.spl')

%EEG1=EEG
EEG=del2map(double(squeeze(mean(ALLEEG(1).data,2))),EEG.chanlocs)



fEEG=eeglab2fieldtrip(ALLEEG(1),'timelockanalysis','none')
cfg.elec         = fEEG.elec ;%structure with electrode positions or filename, see FT_READ_SENS
cfg.trials       = 'all'; 
cfg.feedback     = 'text';%string, 'no', 'text', 'textbar', 'gui' (default = 'text')

%[fff] = ft_senstype(fEEG.elec, 'ant')
%elec = ft_read_sens('D:\WORKINGSET_D\ant64.elp', 'besa_elp') 
[data] = ft_scalpcurrentdensity(cfg, fEEG)

% for the multiple plots also
cfg = [];
cfg.xlim = [0.016:0.002:0.09];
%cfg.ylim = [25 45];
%cfg.zlim = [-1e-27 1e-27];
%cfg.baseline = [-0.5 -0.1];
%cfg.baselinetype = 'absolute';
%cfg.layout = 'CTF151_helmet.mat';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
figure; ft_topoplotTFR(cfg,data);

print(gcf, '-dpdf','-painters','-bestfit',  'MST_rDLPFC_TEP_timeseries_topoplots.pdf')
%% topo/head plot MAX-MIN
   for h=1:size(headminmax,2)
       headminmax2{h}=mean(headminmax{h},3)
   end
   
   figure;hold on; headplot(headminmax2{1}, 'GW64.spl','maplimits',[0 20])% ,'maplimits',[0 0.1]
   figure; headplot(headminmax2{2}, 'GW64.spl','maplimits',[0 20])
   
    figure; hold on; topoplot(headminmax2{1}, EEG.chanlocs,'maplimits',[0 20]);  title([ALLEEG(1).setname], 'interpreter', 'none'); colorbar; hold off;
     figure; hold on; topoplot(headminmax2{2}, EEG.chanlocs,'maplimits',[0 20]); title([ALLEEG(2).setname], 'interpreter', 'none');colorbar; hold off;
%% significance headplot for time window
  clearvars  TEPVec1
datasets=[1 2]
    statwin=[20 40]
    clear TEPVec1
    for i=datasets
    clear statInd
        statInd=ALLEEG(i).times>=statwin(1) & ALLEEG(i).times<=statwin(2);
    %TEPVec1{i}=mean(double(ALLEEG(i).data(:,statInd,:)),1); %mean
    %TEPVec1{i}=sum(abs(double(ALLEEG(i).data(:,statInd,:))),2); %rectfied
    %TEPVec1{i}=trapz(ALLEEG(i).times(statInd),abs(ALLEEG(i).data(:,statInd,:)),2); % AUC
    TEPVec1{i}=sum(abs(ALLEEG(i).data(:,statInd,:)),2)./(ALLEEG(i).srate/1000); % AUC
    end
    TEPVec1=  TEPVec1( find(~cellfun(@isempty,TEPVec1)));
    figure; headplot(mean(TEPVec1{2},3), 'GW64.spl','cbar',[0 120],'maplimits',[0 120]);
    [~, p]=ttest2(TEPVec1{1},TEPVec1{2}, 'Dim', 3);
    p_log=-log10(p)
    figure; headplot(p, 'GW64.spl','cbar',[0 5],'maplimits',[0 0.1])
    figure; headplot(p_log, 'GW64.spl','cbar',[0 4],'maplimits',[0 4]);
    set(gcf, 'Color', 'None')
    figure; topoplot(p_log, 'GW64.spl');colorbar
    %export_fig HEADPLOT_ADHD_CONT_p30.png -q101 -transparent
%% Time series significance topoplot
    clear p
    datasets=[1 2]
    [~, p]=ttest2(ALLEEG(dataset(1)).data,ALLEEG(dataset(2)).data, 'Dim', 3);
    EEG=ALLEEG(1);
    EEGp=EEG;
    EEGp.data(:,:,1)=p;
    EEGp.data(:,:,2)=p;
    EEGp=pop_select(EEGp, 'trial', [1 2]);
    EEGlogp=EEG;
    EEGlogp.data(:,:,1)=-log10(p);
    EEGlogp.data(:,:,2)=-log10(p);
    EEGlogp=pop_select(EEGlogp, 'trial', [1 2]);
    pop_topoplot(EEGp,1, [10:5:50] ,'P significance map',[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0 ,'maplimits',[0 0.1]);
    pop_topoplot(EEGlogp,1, [10:5:50] ,['-log10(P) significance map ' ALLEEG(datasets(1)).condition(8:end) ' compared to ' ALLEEG(datasets(2)).condition(8:end)],[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0,'maplimits',[0 2] );%,'maplimits',[0 0.1]
    %export_fig time_series_heaplots_ADHD_HEALTHY_-log10P.tif -q101 -append
    % print(gcf, '-dpdf','-painters','-bestfit',  'TEP_ECT_PRE_POST.pdf')
    % print(gcf,  'TEP_ECT_PRE_POST.pdf', '-dpdf','-painters',)
    % export_fig time_series_heaplots_ADHD_HEALTHY_-log10P.PNG -q101
    % export_fig time_series_heaplots_Healthy_15ms-40ms-log10P.PNG -q101
    % export_fig time_series_heaplots_Healthy_15ms-40ms-log10P.PDF -q101 -opengl -append
    % averege of time window of the significance/DIFF DATA
%% stat window headplot
    var=EEGi(3)
    statwin=[20 30]
    statInd=var.times>=statwin(1) & var.times<=statwin(2);
    Head=nanmean(nanmean(var.data(:,statInd,:), 3),2);%    TEP(:,:,i)=nanmean(ALLEEG(i).data, 3);
    figure; headplot(Head, 'GW64.spl','electrodes', 'on','labels',2)%,'maplimits',[0 0.1]
    figure; topoplot(Head ,EEG.chanlocs);
    
    % ??
    statwin=[160 250]
    statInd=EEG.times>=statwin(1) & EEG.times<=statwin(2);
    HeadDiff=nanmean(nanmean(EEG.data(:,statInd,:), 3),2);% 
    %figure; headplot(HeadDiff, 'GW64.spl')
    figure;  topoplot(HeadDiff,EEG.chanlocs)
       pop_topoplot(EEG,1, [10:5:250] ,'significance map',[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0); %,'maplimits',[-5 5]
    figure; headplot(HeadDiff, 'GW64.spl','electrodes', 'off') %'maplimits',maplim
    
   % export_fig Headplot_ADHD_HEALTHY_DIFF_P30_cord_165_12.PNG -q101
%% Cohen's D time windows mean topoplot 
% Cohen's D time windows mean topoplot 
alpha=0.05
% for a grand file
clear datasets EEG1 EEG2 timeWin timeind eeg1 eeg2 pool_std t_struct pool_var_numertor pool_var_dom pool_var pool_std
datasets=[1 2]
EEG1=ALLEEG(datasets(1));
EEG2=ALLEEG(datasets(2));
% export_fig cohens_D_topoplot_p60_ADHD_compared_HEALTHY_15ms_compared_to_40ms_cohens_D.PNG -q101
timeWin=[20 40]
maplim=[0 1] %topoplot maplimits
timeind=ALLEEG(1).times>=timeWin(1) & ALLEEG(1).times<=timeWin(2);
eeg1 = squeeze(mean(EEG1.data(:,timeind,:),2)); % egg1,2: (electrodes, subjects)
eeg2 = squeeze(mean(EEG2.data(:,timeind,:),2));
%for head matrix (head / headminmax computed be TEPstats)
%eeg1=squeeze(head{1})
%eeg2=squeeze(head{2});
%pool_std = std([eeg1 eeg2],1,2);
pool_var_numertor = (std(eeg1, 0, 2).^2*(size(eeg1,2)-1))+(std(eeg2, 0, 2).^2*(size(eeg2,2)-1));
pool_var_dom = size(eeg1,2) + size(eeg2,2)-2;
pool_var = pool_var_numertor./pool_var_dom;
pool_std = sqrt(pool_var); 
 % Cohens_D
[t_struct.h,t_struct.p] = ttest2(eeg1,eeg2,'alpha',alpha, 'Dim', 2);
[t_struct.effect_size] = (mean(eeg1, 2) - mean(eeg2, 2))./pool_std;
%load a eeglab dataset that contains EEG.chanlocs
figure; topoplot(t_struct.effect_size, EEG.chanlocs,'maplimits',maplim,'conv','on','style','map','shading','interp'); 
colorbar; title([{'Cohen''s D' , ALLEEG(datasets(1)).setname , ALLEEG(datasets(1)).setname}]);
set(get(gca,'title'),'Position',[0 0.45 1.00011]);
% export_fig ADHD_plus_HEALTHY_15ms_compared_to_40ms_cohens_D.PNG -q101
%% Cohen's D Cohens_D time series
alpha=0.05
group='withi' % dependet groups='within' independent=''
% for a grand file
clear datasets EEG1 EEG2 timeWin timeind eeg1 eeg2 pool_std t_struct pool_var_numertor pool_var_dom pool_var pool_std
datasets=[2 1]
EEG1=ALLEEG(datasets(1));
EEG2=ALLEEG(datasets(2));

eeg1=EEG1.data;%(:,:,[1:45]);
eeg2=EEG2.data;%(:,:,[1:45]);
% eeg1=cat(3,ALLEEG(1).data,ALLEEG(2).data); % all datasets - not grouped
% eeg2=cat(3,ALLEEG(3).data,ALLEEG(4).data); % all datasets - not grouped
pool_var_dom = size(eeg1,3) + size(eeg2,3)-2;

for i=1:size(eeg1,2)
    pool_var_numertor = (std(squeeze(eeg1(:,i,:)), 0, 2).^2*(size(squeeze(eeg1(:,i,:)),2)-1))+(std(squeeze(eeg2(:,i,:)), 0, 2).^2*(size(squeeze(eeg2(:,i,:)),2)-1));
    pool_var = pool_var_numertor./pool_var_dom;
    pool_stdT = sqrt(pool_var);
    if strcmpi(group,'within') && size(eeg1,3)==size(eeg2,3)
        r = diag(corr(squeeze(eeg1(:,i,:))',squeeze(eeg2(:,i,:))'));
        pool_std(:,i) = pool_stdT.*sqrt((1-r).*2);
    else
        pool_std(:,i) = pool_stdT;
    end
end
[t_struct.h,t_struct.p] = ttest2(eeg1,eeg2,'alpha',alpha, 'Dim', 3);
[t_struct.effect_size] = (mean(eeg1, 3) - mean(eeg2, 3))./pool_std;

coh=ALLEEG(datasets(1)); coh.data=[];
coh.data=t_struct.effect_size;
coh.times=ALLEEG(datasets(1)).times;
coh.xmin=(coh.times(1)/1000); coh.xmax=(coh.times(1,end)/1000);
coh.srate=ceil(size(coh.data,2)/(coh.xmax-coh.xmin))-1; % (xmax-xmin)*srate+1 = number of frames
coh.chanlocs=EEG1.chanlocs;
coh.data=[];
coh.data(:,:,1)=t_struct.effect_size; coh.data(:,:,2)=t_struct.effect_size;
coh=pop_select(coh, 'trial', [1 2]);
pop_topoplot(coh,1, [10:5:150] ,'cohens D',[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0); %,'maplimits',[-5 5]
%pop_topoplot(coh,1, [45:5:250] ,'cohens D',[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0); %,'maplimits',[-5 5]
% export_fig ADHD_plus_HEALTHY_15ms_compared_to_40ms_cohens_D.PNG -q101
%% supfunsim -  reconstruction of sources of brain activity and reconstruction of EEG signa
addpath(genpath(['D:\Google_drive\MATLAB\supFunSim']));
parameters = EEGParameters().generate();
filters = [ "LCMV", "MMSE", "ZF", "RANDN", "ZEROS", "EIG-LCMV", ...
            "sMVP_MSE", "sMVP_R", "sMVP_N", "sMVP_NL_MSE", ...
            "sMVP_NL_R", "sMVP_NL_N", "NL" ];
reconstruction = EEGReconstruction();
reconstruction = reconstruction.init();
for np = 1:length(parameters)
    parameter = parameters(np);
    reconstruction = reconstruction.setparameters(parameter);
    reconstruction = reconstruction.setsignals();
    reconstruction = reconstruction.setleadfields();
    reconstruction = reconstruction.setpreparations();
    reconstruction = reconstruction.setfilters(filters);
    reconstruction.save();
end
reconstruction = reconstruction.printaverageresults();
eegplot = EEGPlotting(reconstruction)
eegplot.plotgausswave()
eegplot.plotMVARzeroingmatrix()
eegplot.plotMVARmodelcoefficientmatrix2()
eegplot.plotPDC2()
%% just regular headplots
    maplimits=[-0.5 1]
     figure; hold on; headplot(TEPVecstat(:,:,2), 'GW64.spl','electrodes', 'off')%,'maplimits',maplimits)
    scatter3(EEG.chanlocs(2).X, EEG.chanlocs(2).Y, EEG.chanlocs(2).Z, 'k','filled')
    
    figure; headplot(TEPVecstat(:,:,1), 'GW64.spl','electrodes', 'on', 'transform',[-0.35579 -6.3369 12.3705 0.053324 0.018746 -1.5526 1.0637 0.98772 0.93269])
    figure; headplot(TEPVecstat(:,:,6), 'GW64.spl','electrodes', 'off') %'maplimits',maplimits)
    
    figure; headplot(TEPVecstat(:,:,6), 'GW64.spl','electrodes', 'on','maplimits',maplimits)
    figure; headplot(TEPVecstat(:,:,6), 'GW64.spl','electrodes', 'on','maplimits',maplimits)
%% old version of p-value headplot
    %file=2 %/2
    statInd=EEGp.times>=statwin(1) & EEGp<=statwin(2);
    HeadTtest1(:,:)=nanmean(ALLEEG(1).data(:,statInd,:), 2);
    HeadTtest2(:,:)=nanmean(ALLEEG(2).data(:,statInd,:), 2);
    
    for e=1:64
        [H(e), Headt(e), ~, stats]=ttest2(log10(abs(HeadTtest1(e,:))),log10(abs(HeadTtest2(e,:))))
    end
    
    figure; headplot(Headt, 'GW64.spl','electrodes', 'off') %'maplimits',maplimits)
%% topoplot series??
    dat=6
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',dat,'study',0);
    pop_topoplot(EEG,dat, [15:10:250] ,ALLEEG(dat).setname,[4 6] ,0,'electrodes','on','maplimits',[-4 4]); % [statwin(1):10:statwin(2)]
    export_fig topoplot_series.tif -q101 -append -opengl
    %close all
%% Wavelet per subject    
% electrodes={'f5' 'f3' 'fc5' 'fc3' 'fcz' 'cz' 'fcz'}; % {'ALL'}
% electrodes={'f5'}; % {'ALL'}
% %elecmethod='no';
% T1=-2000;
% T2=2000;
% %erspmax= 2;
% %itcmax= 0.3;
% %padratio= 16;
% elocFile= [];
% plotgraph=  'off' ;
% savefile= 'on';
% %baseWin= [-1500 -600];
% % method= 'div';
% % trialbase= 'on';
% freqs=[2 56];
%TMS_NewWavelet2(electrodes, elecmethod, T1, T2, freqs, erspmax, itcmax, padratio, elocFile, plotgraph, savefile, baseWin, method, trialbase)
%TMS_NewWavelet2(electrodes, T1, T2, freqs, elocFile, plotgraph, savefile)
%TMS_NewWavelet2

TMS_NewWavelet3
%% Wavelet GRAND files
% remember to check baseline correction
grandName='PRE';
bc=0;
%wavelet(j,1).tfdata_blc=Z_baseLineCorrect(wavelet(j,1).tfdata, baseWin, wavelet(j,1).times, method); 
TMS_WaveGrandAv(grandName)

% draw_specto2
%% JUNK
    
    
    
    a = [150 250];
    maplimits=[0 1]
    EEG = pop_headplot(EEG, 0, a , ['Components of dataset: ' EEG.setname], [], 'setup',{['GW64.spl'] 'meshfile' 'mheadnew.mat' 'transform' [-0.35579 -6.3369 12.3705 0.053324 0.018746 -1.5526 1.0637 0.98772 0.93269] }, 'electrodes', 'off', 'view', [180 0]);
    
    
    
    maplimits=[0 1]
    figure; topoplot(ERPVecstat(:,:,1), EEG.chanlocs, 'electrodes','on', 'style', 'map' ,'shading','interp','headrad',0)%,'maplimits',maplimits); %  'maplimits',maplimits
    EEG.data_temp=EEG.data
    EEG.data=ERPVecstat(:,:,1)
    figure; headplot(EEG, EEG.chanlocs, 'electrodes','on')%,'maplimits',maplimits); %  'maplimits',maplimits
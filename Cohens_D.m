% Cohens_D
%% Cohen's D (time series topoplots / mean time win)
% computing the p value and the cohen's D effect size - FOR UNPAIRED sample.
% open Cohens_D
alpha=0.05
% for a grand file
clear datasets EEG1 EEG2 timeWin timeind eeg1 eeg2 pool_std t_struct pool_var_numertor pool_var_dom pool_var pool_std
datasets=[1 3]
EEG1=ALLEEG(datasets(1));
EEG2=ALLEEG(datasets(2));
timeind=ALLEEG(1).times>=timeWin(1) & ALLEEG(1).times<=timeWin(2);
%% for a grand file
EEG1=ALLEEG(1);
EEG2=ALLEEG(2);
eeg1 = squeeze(mean(EEG1.data(:,timeind,:),2)); % egg1,2: (electrodes, subjects)
eeg2 = squeeze(mean(EEG2.data(:,timeind,:),2));
%% OR for head matrix (head / headminmax computed be TEPstats)
eeg1=squeeze(head{1});
eeg2=squeeze(head{2});
%% Cohen's D time windows mean topoplot 
% Cohen's D time windows mean topoplot 
timeWin=[15 300]
timeind=ALLEEG(1).times>=timeWin(1) & ALLEEG(1).times<=timeWin(2);
eeg1 = squeeze(mean(EEG1.data(:,timeind,:),2)); % egg1,2: (electrodes, subjects)
eeg2 = squeeze(mean(EEG2.data(:,timeind,:),2));
%for head matrix (head / headminmax computed be TEPstats)
eeg1=squeeze(head{1});
eeg2=squeeze(head{2});
%pool_std = std([eeg1 eeg2],1,2);
pool_var_numertor = (std(eeg1, 0, 2).^2*(size(eeg1,2)-1))+(std(eeg2, 0, 2).^2*(size(eeg2,2)-1));
pool_var_dom = size(eeg1,2) + size(eeg2,2)-2;
pool_var = pool_var_numertor./pool_var_dom;
pool_std = sqrt(pool_var); 
 % Cohens_D
[t_struct.h,t_struct.p] = ttest2(eeg1,eeg2,'alpha',alpha, 'Dim', 2);
[t_struct.effect_size] = (mean(eeg1, 3) - mean(eeg2, 3))./pool_std;
%load a eeglab dataset that contains EEG.chanlocs
figure; topoplot(t_struct.effect_size, EEG.chanlocs,'maplimits',[0 1.5],'style','map','shading','interp'); colorbar; title('Cohen''s D')
%% Cohen's D Cohens_D time series
eeg1=EEG1.data;
eeg2=EEG2.data;
% eeg1=cat(3,ALLEEG(1).data,ALLEEG(2).data); % all datasets - not grouped
% eeg2=cat(3,ALLEEG(3).data,ALLEEG(4).data); % all datasets - not grouped

pool_var_dom = size(eeg1,3) + size(eeg2,3)-2;

for i=1:size(eeg1,2)    
r = diag(corr(squeeze(eeg1(:,i,:))',squeeze(eeg2(:,i,:))'));
pool_var_numertor = (std(squeeze(eeg1(:,i,:)), 0, 2).^2*(size(squeeze(eeg1(:,i,:)),2)-1))+(std(squeeze(eeg2(:,i,:)), 0, 2).^2*(size(squeeze(eeg2(:,i,:)),2)-1));
pool_var = pool_var_numertor./pool_var_dom;
pool_stdT = sqrt(pool_var); 
pool_std(:,i) = pool_stdT.*sqrt((1-r).*2);
end
[t_struct.h,t_struct.p] = ttest2(eeg1,eeg2,'alpha',alpha, 'Dim', 3);
[t_struct.effect_size] = (mean(eeg1, 3) - mean(eeg2, 3))./pool_std;

cohens=load('D:\Itay-Google-Drive\MATLAB\LAB_MatlabScripts\EEGlab_empty_stracture.mat');
coh=cohens.EEG;
coh.data=t_struct.effect_size;
coh.times=ALLEEG(datasets(1)).times;
coh.xmin=(coh.times(1)/1000); coh.xmax=(coh.times(1,end)/1000);
coh.srate=ceil(size(coh.data,2)/(coh.xmax-coh.xmin))-1; % (xmax-xmin)*srate+1 = number of frames
coh.data=[];
coh.data(:,:,1)=t_struct.effect_size; coh.data(:,:,2)=t_struct.effect_size;
coh=pop_select(coh, 'trial', [1 2]);
pop_topoplot(coh,1, [10:5:250] ,'cohens D',[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0); %,'maplimits',[-5 5]
pop_topoplot(coh,1, [45:5:250] ,'cohens D',[5 10] ,0,'electrodes','on', 'style', 'map', 'headrad', 0); %,'maplimits',[-5 5]

%% plot
%load a eeglab dataset that contains EEG.chanlocs
figure; topoplot(t_struct.effect_size, EEG.chanlocs,'maplimits',[0 1.5],'style','map','shading','interp'); colorbar; title('Cohen''s D')
% export_fig ADHD_HEALTHY_MIN_MAX_cohens_D.PNG -q101
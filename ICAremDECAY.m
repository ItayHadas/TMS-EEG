function [EEG] = ICAremDECAY(EEG,ROI)


% ROI = {'Fpz' 'Fp2' 'AF8' 'F8' 'FT8' 'AFz' 'AF4' 'F6' 'FC6' 'Fz' 'F2' 'F4' 'FC4' 'FC2' };
ROIThresh=0.65;
%ROI = {'F3' 'FP1' 'fpz' 'fp2' 'AF3' 'F11' 'FT11' 'F12' 'FT12' 'F1' 'FT7' 'AF5' 'F7' 'FC7' 'FC5' 'F5' 'FC3' 'c3' 'c5' 't7' 'af7' 'fpz' 'afz' 'fz' 'fcz'};
tmsMuscleThresh=5; % best threshold is 7 if the the decay is longer try to reduce the number
maxcomps=15 ;%how many ICA components to check
remocomp=5 ;% how many ICA comp. to remove.

%% loops through files -omited
% for i=1: size(fileNames,1)
%     if size(fileNames, 1)==1
%         fileName=fileNames{i,1}';
%     else
%         fileName=fileNames{i,1};
%     end
%      [ num2str(i) ' out of ' num2str(size(fileNames,1)) '   -   ' fileName ]
%     EEG = pop_loadset( ['D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\ICA1\WEL003_SPD_BL_Epo_Dec_interp_ICA1.set']);
%D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\ICA1\WEL003_SPD_BL_Epo_Dec_interp_ICA1.set
%     EEG = eeg_checkset( EEG );
%     %% report table
%     out = table;
%     out.Subject = {[EEG.subject]};
%     out.filename = {EEG.filename};
%     setname= EEG.setname;
%     datfile=EEG.datfile; 

%%
warning 'off'
%Ranks components, outputs variance
[EEG, ~] = tesa_sortcomps(EEG);
comps = 1:maxcomps; % 1:size(EEG.icaweights,1);
eNum= sort(elecName(EEG, ROI));
% TMS_pop_selectcomps(EEG, [1:35] );
%Calculates time course and variance
compTimeCourse = arrayfun(@(x)eeg_getdatact(EEG, 'component', [x], 'projchan', [eNum]),1:size(EEG.icawinv,2),'UniformOutput', false);
% pop_selectcomps(EEG, [1:28] );
%figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo)
for compNum=comps % compNum=1
    tempCompZ = zscore(EEG.icawinv(:,compNum));
    ActivRatio = mean(abs(tempCompZ(eNum,:)));
    suspects(compNum)= [ActivRatio >= ROIThresh];
end
clear -var tempCompZ ActivRatio compNum

compRem=zeros(size(EEG.icaweights,1),1);
CompWave=[];
comps=comps(suspects);
%if size(comps,2)>1 ; comps=comps(:,1:remocomp); end
% tmsMuscleThresh=2
for compNum=comps  % (1:remocomp) %compNum=comps(suspects(1))
    temp = compTimeCourse{1,compNum};
    [~, elec]=max(abs(EEG.icawinv(eNum,compNum))); 
    temp = squeeze(temp(elec,:,:));
    comp = abs(mean(temp,2));
    [~, t0] = min(abs(EEG.times-(0)));
    [~, t50] = min(abs(EEG.times-(50)));
    % [~, tcut] = min(abs(EEG.times-(EEG.tmscut(1).cutTimesTMS(2))));
    winScore = max(comp(t0:t50,:));
    %nonwinScore=  mean([comp(find(EEG.times>-350 & EEG.times<-50));  comp(find(EEG.times>50 & EEG.times<350))]);
    nonwinScore=  mean([comp([EEG.times>-350 & EEG.times<-50]);  comp([EEG.times>50 & EEG.times<350])]); %changed 20200212 - should work
    % mean([comp(find(EEG.times<-50));  comp(find(EEG.times>300 ))])
    tmsMuscleRatio = winScore./nonwinScore;
    
    if tmsMuscleRatio>tmsMuscleThresh & winScore>25
        compRem(compNum)=1;
        CompWave(end+1,:)=comp';
        %a(end+1,:)=zscore(diff(comp'));
    end
    clear -var temp comp winScore nonwinScore tmsMuscleRatio elec
end
%    [find(compRem,remocomp)']
% pop_plotdata(EEG, 0, [1:size(EEG.icaweights,1)] , [1:size(EEG.data,3)], EEG.setname, 0, 1, [0 0]);
% pop_selectcomps(EEG, [1:28] );

if sum(compRem)>0
    setname1=EEG.setname;
    filename1=EEG.filename;
    datfile1=EEG.datfile;
    EEG=pop_subcomp( EEG, [find(compRem,remocomp)'], 0);
    EEG.setname=setname1;
    EEG.filename=filename1;
    EEG.datfile=datfile1;
    %% printing removed component
    %if figs==1
    %figure; hold on; plot(EEG.times,CompWave([1],:),'-b'); %removed component wave plotting.
    %if size(CompWave,1)>1 ; plot(EEG.times,CompWave([2],:),'-r'); end ;
    % figure; hold on; plot(EEG.times,CompWave(3,:))
    %legend(['subject:' num2str(EEG.subject) '  condition:' num2str(EEG.condition) '  session:' num2str(EEG.session) '  comp: ' num2str(find(compRem,remocomp)')] ); hold off;
    %         end
    %         % report table
    %         out.number_of_muscle_exponent_components_detected = {num2str(sum(compRem))};
    %         out.number_of_muscle_exponent_components_rejected = {num2str(size(find(compRem,remocomp),1))};
    %         out.muscle_exponent_components_detected = {num2str(find(compRem))}; % {num2str(selectedtrials)};
    %         out.muscle_exponent_components_rejected = {num2str(find(compRem,1))}; %
    %         EEG.filename=fileName;
    %         EEG.setname=setname;
    %        EEG.datfile=datfile;
    %%
    EEG.Removed_ICA_Comps.EXPwaveform=CompWave(1:size(find(compRem,remocomp),1),:);
    disp('Exponent components removed:')
    EEG.Removed_ICA_Comps.comps_num={(find(compRem,remocomp))} %

    EEG=Z_append(EEG,'_EXP');
else
    %% report table
    %         out.number_of_muscle_exponent_components_detected = {'non detected'};
    %         out.number_of_muscle_exponent_components_rejected = {'non detected'};
    %         out.muscle_exponent_components_detected = {'non detected'}; % {num2str(selectedtrials)};
    %         out.muscle_exponent_components_rejected = {'non detected'}; %  %
    %%
    EEG=Z_append(EEG,'_NOexp');
    
end
%EEG.filename=cellstr(out.filename);
%EEG.setname=setname;
%EEG=Z_append(EEG,' auto_exp_prune');

%EEG.Removed_ICA_Comps.High_Amp_Exp={num2str(find(compRem,1))};
EEG = eeg_checkset( EEG );
%EEG = pop_saveset( EEG, [pathName 'ICA1\' EEG.filename]);
end
%% table - ommited
%output(i,:)=out;

%clear -var EEG compTimeCourse comps compRem suspects compNum out eNum CompWave
% end
%
% %% report table
% output=sortrows(output,1)
% warning 'on'
% writetable(output(:,1:6),strjoin([pathName 'output_ICA_exp_comp_subjects-' output.Subject(1,:) '-' output.Subject(size(output.Subject,1),:) '.xlsx'])) ;
% %%
% clearvars -except output;
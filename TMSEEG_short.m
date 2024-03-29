%% preprocess pipe
% pathi='D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\';
% [fileNames, pathName]=Z_getSetsFileNames('set',pathi);
clear all
%[fileNames, pathName]=Z_getSetsFileNames('cdt');
outdir='A:\WorkingSet\suicidality_TEP\' %pathName
Path = dir('\\ad.ucsd.edu\ahs\apps\INTERPSYC\DATA\Suicidality_801566\Neurophysiology_Data\**\*SPD_*.set');
fileNames={Path.name}';

for i=1: size(fileNames,1) 
    disp(['***************  dataset ' num2str(i) '/' num2str(size(fileNames,1)) '  **************']);
    fileName=fileNames{i};
    pathName=[Path(i).folder '\'];
    EEG = pop_loadset( [pathName fileName]);
    trim=[-2 15]; %DoublePulseINT=0 % trim=[-2 5];
    EEG = Triminterp(EEG,trim);
    base=[-200   -20];
    EEG  = pop_rmbase( EEG, base );
    EEG.data=double(EEG.data);
    EEG = pop_runica(EEG, 'icatype', 'picard','m',12 , 'maxiter',700); %,'mode','standard';
    EEG = Z_append(EEG,'_ICA1');
    mkdir([outdir 'ICA1\'])
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'ICA1\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
    EEG = ICAremDECAY(EEG);
    mkdir([pathName '\removed_comps_ICA1'])
    minfreq=58; maxfreq=62;
    EEG = FFT_comp(EEG, minfreq, maxfreq);
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'removed_comps_ICA1\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
    low=0.2;    high=80;
    EEG.data = FilterTEP(EEG.data, low, high,EEG.srate); 
    EEG.data=double(EEG.data);
    for elecc=1:size(EEG.data,1)
        for epoo=1:size(EEG.data,3)
            EEG.data(elecc,:,epoo)=detrend(double(EEG.data(elecc,:,epoo)),1);
        end
    end
    EEG = eeg_checkset( EEG );  
    %EEG = pop_reref(EEG, []);
    base=[-800  -200];
    EEG  = pop_rmbase( EEG, base );
    EEG.data=double(EEG.data);
    %EEG = pop_runica(EEG,'icatype','fastica','approach','symm','g','tanh','stabilization','on');
    EEG = pop_runica(EEG, 'icatype', 'picard','m', 12 , 'maxiter',700); %,'mode','standard';
    EEG = Z_append(EEG,'_ICA2');
    mkdir([pathName '\ICA2'])
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'ICA2\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
    clearvars -except chanlocs i fileNames pathName outdir Path
end

%% remove comps
clear all
%pathi='D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\';
pathi='A:\WorkingSet\suicidality_TEP\ICA2\';
[fileNames, pathName]=Z_getSetsFileNames('set',pathi);
i=0
%for i=1: size(fileNames,1) 
%%
i=i+1
    fileName=fileNames{i};
    EEG = pop_loadset( [pathName fileName]);
    ALLEEG(1)=EEG;
    CURRENTSET=1;
    %  pop_eegplot( EEG, 1, 1, 1);
    [EEG, ~] = sortcomps(EEG, [0 500]);
    EEG = pop_iclabel(EEG, 'default');
    pop_selectcomps(EEG, [1:28] );
    %pop_viewprops( EEG, 0, [1:28], {'freqrange', [2 60]}, {}, 1, '' )
    
%%
if (contains(EEG.setname,'_RS'))
channels= elecName(EEG,{'f4'}) %41  27 47 20
elseif (contains(EEG.setname,'_LS'))
    channels= elecName(EEG,{'f3'}) %41  27 47 20
end
%components=[1     2     3     5    13];
%components= [ 1 2 ] % components= RejComp'
components= [find(EEG.reject.gcompreject)] %3     4     5     6     7     8     9    11    35
comp='off';
comprem(EEG,channels,components,comp)
    %RejComp=sort
    % ([find(any(EEG.etc.ic_classification.ICLabel.classifications(:,4:6)>0.9,2));...
    %    find(any(EEG.etc.ic_classification.ICLabel.classifications(:,2:3)>0.86,2))]);
    %EEG = pop_subcomp( EEG,RejComp , 0);
    %EEG = Z_append(EEG,'_remcomp');
%%
EEG = pop_subcomp( EEG, components, 0);
EEG = Z_append(EEG,'_final1');
mkdir([pathName '\final1'])
EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([pathName 'final1\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
close all
disp(['*************** Saved dataset ' num2str(i) '/' num2str(size(fileNames,1)) '**************']);

%% GRAND Average 

%load('D:\MATLAB\LAB_MatlabScripts\Chanlocs\chanlocs66_flexnet_compumedics.mat'); 
%pathi='D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\ICA2\final1\'
pathi='A:\WorkingSet\suicidality_TEP\ICA2\final1\'
chan_interp='on'
chanlocs=chanlocs66;
STDcalc=0
Z2_grand_average('SUI_POST_LSP',[1:6],pathi, chan_interp,chanlocs,STDcalc) 

% figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:))); ylim([-40 40]); xlim([-300 300])

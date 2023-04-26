%%
clear all
[fileNames, pathName]=Z_getSetsFileNames('set');

for i=1: size(fileNames,1) 
    fileName=fileNames{i};
    EEG = pop_loadset( [pathName fileName]);
    base=[-200   -20];
    EEG  = pop_rmbase( EEG, base );
    EEG = pop_runica(EEG, 'icatype', 'picard', 'maxiter',500); %,'mode','standard';
    EEG = Z_append(EEG,'_ICA1');
    mkdir([pathName '\ICA1'])
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([pathName 'ICA1\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
    EEG = ICAremDECAY(EEG);
    mkdir([pathName '\removed_comps_ICA1'])
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([pathName 'removed_comps_ICA1\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
    low=0.2;    high=80;
    EEG.data = FilterTEP(EEG.data, low, high,EEG.srate); 
    EEG = pop_reref(EEG, []);
    EEG.data=double(EEG.data);
    for elecc=1:size(EEG.data,1)
        for epoo=1:size(EEG.data,3)
            EEG.data(elecc,:,epoo)=detrend(double(EEG.data(elecc,:,epoo)),1);
        end
    end
    EEG = eeg_checkset( EEG );  
    base=[-800   -200];
    EEG  = pop_rmbase( EEG, base );
    EEG = pop_runica(EEG, 'icatype', 'picard', 'maxiter',500);
    EEG = Z_append(EEG,'_ICA2');
    mkdir([pathName '\ICA2'])
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([pathName 'ICA2\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
    clearvars -except chanlocs i fileNames pathName
end 

%% remove comps
clear all
[fileNames, pathName]=Z_getSetsFileNames('set');

for i=1: size(fileNames,1) 
    fileName=fileNames{i};
    EEG = pop_loadset( [pathName fileName]);
    [EEG, ~] = tesa_sortcomps(EEG);
    EEG = pop_iclabel(EEG, 'default');
    RejComp=sort([find(any(EEG.etc.ic_classification.ICLabel.classifications(:,4:6)>0.9,2));...
        find(any(EEG.etc.ic_classification.ICLabel.classifications(:,2:3)>0.86,2))]);
    EEG = pop_subcomp( EEG,RejComp , 0);
    EEG = Z_append(EEG,'_remcomp');
    mkdir([pathName '\remcomp'])
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([pathName 'remcomp\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
end

%% GRAND Average 

%load('D:\MATLAB\LAB_MatlabScripts\Chanlocs\chanlocs66_flexnet_compumedics.mat'); 

chan_interp='on'
chanlocs=chanlocs66;
STDcalc=0
Z2_grand_average('WL_SP_all',[1:6],chan_interp,chanlocs,STDcalc) 

figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:))); ylim([-40 40]); xlim([-300 300])

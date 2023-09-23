clear all
addpath('D:\GITs\AARATEPPipeline');
outdir='A:\WorkingSet\suicidality_TEP\' %pathName
Path = dir('\\ad.ucsd.edu\ahs\apps\INTERPSYC\DATA\Suicidality_801566\Neurophysiology_Data\**\*SPD_*.cdt');
fileNames={Path.name}';
chanlocs62=load('D:\MATLAB\LAB_MatlabScripts\Chanlocs\chanlocs66_flexnet_compumedics.mat');
%  chanlocs = readlocs( 'd:\\Google_Drive\\MATLAB\\EEGLAB\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
%  load('D:\Google_Drive\MATLAB\LAB_MatlabScripts\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat');%ST-THETA-BURST
for i=1: size(fileNames,1) %%%%%%%%%%%%%%%%% PIPELINE LOOP
    fileName=fileNames{i};
    pathName=[Path(i).folder '\'];
    if strcmpi(fileName(:,end-2:end),'set')
        EEG = pop_loadset( [pathName fileName]);
    elseif strcmpi(fileName(:,end-2:end),'cnt')
        EEG = pop_loadcnt([pathName fileName] , 'dataformat', 'int32');
    elseif strcmpi(fileName(:,end-2:end),'cdt')
        EEG = loadcurry([pathName fileName], 'CurryLocations', 'True'); %'CurryLocations', 'False'
    elseif strcmpi(fileName(:,end-3:end),'vhdr')
        EEG  = pop_loadbv(pathName , fileName );
    end
    disp([ num2str(i) ' out of ' num2str(size(fileNames,1)) '   -   ' fileName ])
    EEG.filename=fileName;
    EEG.subject=EEG.filename(1:6);
    if strcmpi(fileName(:,end-3:end),'vhdr')
        EEG.setname=EEG.filename(1:end-5);
    else
        EEG.setname=EEG.filename(1:end-4);
    end
    % pop_eegplot( EEG, 1, 1, 1);
    if contains(EEG.setname,'BL','IgnoreCase',true)
        EEG.condition='pre';
    elseif contains(EEG.setname,'post','IgnoreCase',true)
        EEG.condition='post';
    end
    EEG.session='1';
    EEG.data=double(EEG.data);
    EEG = eeg_checkset( EEG );
    disp(['Dataset is loaded ' num2str(i) '/' num2str(size(fileNames,1))])
    disp(['remove un-needed channels  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    remove_channels= {'lz' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo)
    chanlocs62=EEG.chanlocs;
    EEG = pop_select( EEG,'nochannel',remove_channels);
    EEG = eeg_checkset( EEG );
    for i=1:size(EEG.event,2)
        EEG.event(i).type = num2str(EEG.event(i).type);
    end
    EEG = c_TMSEEG_Preprocess_AARATEPPipeline(EEG,...
        'epochTimespan', [-0.5 0.5],...
        'outputDir', outdir,...
        'outputFilePrefix', 'Clean',...
        'doDecayRemovalPerTrial',true);  
matlab -nosplash -nodesktop

workspace 
clear all
if (ispc)
    sep='\';
    not_sep='/';
    rep_space = ' ';
    GITS='D:\GITs\';
    Path = dir('\\ad.ucsd.edu\ahs\apps\INTERPSYC\DATA\Wellcome_Leap_802232\Neurophysiology_Data\**\*SPD_*.cdt');
    outdir='A:\WorkingSet\WellcomeLeap_TEP';
elseif (ismac || isunix)
    sep='/';
    not_sep='\';
    rep_space = '\ ';
    GITS='/media/ipp/DATA/GITs';
    Path = dir('/mnt/INTERPSYC/DATA/Wellcome_Leap_802232/Neurophysiology_Data/**/*SPD_*.cdt');
    outdir='/media/ipp/DATA/EEG_DATA/WellcomeLeap_TEP';
    eeglabdir='/media/ipp/DATA/Documents/MATLAB/eeglab2023.1'; 
end

addpath([GITS sep 'TMS-EEG']); 
addpath([GITS sep 'TMS-EEG' sep 'picard']); 
addpath([GITS sep 'TESA']);
addpath(genpath([GITS sep 'AARATEPPipeline' sep 'Common' sep 'ThirdParty' sep 'amica']));
gitdir=[GITS sep 'AARATEPPipeline'];
addpath(gitdir);
addpath([GITS sep 'AARATEPPipeline' sep 'Common']);
addpath([GITS sep 'AARATEPPipeline' sep 'Common' sep 'EEGAnalysisCode']);
addpath(eeglabdir)

fileNames={Path.name}';
%chanlocs=load('D:\MATLAB\LAB_MatlabScripts\Chanlocs\chanlocs66_flexnet_compumedics.mat');
eeglab nogui
for i= 1: size(fileNames,1) %%%%%%%%%%%%%%%%% PIPELINE LOOP
    fileName=fileNames{i};
    pathName=[Path(i).folder sep];
    disp([ 'Loading ' num2str(i) ' out of ' num2str(size(fileNames,1)) '   -   ' fileName ])
    if strcmpi(fileName(:,end-2:end),'set')
        EEG = pop_loadset( [pathName fileName]);
    elseif strcmpi(fileName(:,end-2:end),'cnt')
        EEG = pop_loadcnt([pathName fileName] , 'dataformat', 'int32');
    elseif strcmpi(fileName(:,end-2:end),'cdt')
        EEG = loadcurry([pathName fileName], 'CurryLocations', 'True');
    elseif strcmpi(fileName(:,end-3:end),'vhdr')
        EEG  = pop_loadbv(pathName , fileName );
    end
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
    remove_channels= { 'F11' 'F12' 'FT11' 'FT12' 'CB1' 'CB2' 'lz' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo)
    %chanlocs62=EEG.chanlocs;
    EEG = pop_select( EEG,'nochannel',remove_channels);
    evname= '128'; 
    newtrigs=1; 
    method=1;
    DoublePulseINT=0; % paired-pulse interval, Zero (0) in case of single pulse
    [impotant, ~]=elecName(EEG,{'f5' 'f3' 'f1' 'fc5' 'fc3' 'fc1' 'fcz' 'fc2' 'c1' 'cz' 'c2' 'c4' 'c6' 'cp5' 'cp3'});
    if (contains(EEG.setname,'_RS'))
        [impotant, ~]=elecName(EEG,{'f6' 'f4' 'f2' 'fc6' 'fc4' 'fc1' 'fcz' 'fc2' 'c1' 'cz' 'c2' 'c4' 'c6' 'cp5' 'cp3'});
    end
    EEG = IDpulse(EEG,method,DoublePulseINT,impotant,newtrigs,evname);
    EEG = eeg_checkset( EEG );
    EEG_mat = c_TMSEEG_Preprocess_AARATEPPipeline(EEG,...
        'epochTimespan', [-1 1],...
        'outputDir', outdir,...
        'pulseEvent',evname,...
        'outputFilePrefix',[EEG.setname '_Clean'],...
        'ICAType','fastica',...
        'doDecayRemovalPerTrial',true);   %'ICAType','amica',...
    EEG_mat.setname=[EEG.setname '_AARATEPPipeline'];
    EEG_mat.filename=[EEG.setname '_AARATEPPipeline.set']
    EEG_mat.filepath=outdir;
    EEG_mat.datfile='';
    EEG_mat = pop_saveset( EEG_mat, 'filename',EEG_mat.filename,'filepath',EEG_mat.filepath,'check', 'on','savemode','onefile','version','7.3');
    clear EEG_mat EEG
end
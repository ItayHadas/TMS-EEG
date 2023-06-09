clear all

[fileNames, pathName]=Z_getSetsFileNames('cdt');
% [fileNames, pathName]=Z_getSetsFileNames('set');aa
%[fileNames, pathName]=Z_getSetsFileNames('cdt');
outdir='A:\WorkingSet\suicidality_TEP\' %pathName
Path = dir('\\ad.ucsd.edu\ahs\apps\INTERPSYC\DATA\Suicidality_801566\Neurophysiology_Data\**\*SPD_*.cdt');
fileNames={Path.name}';
% [fileNames, pathName]=Z_getSetsFileNames('set');
%  pathName='D:\WORDKINGset_D\SARA_theta_jittered_data\'
%  chanlocs62=load('D:\MATLAB\LAB_MatlabScripts\Chanlocs\chanlocs66_flexnet_compumedics.mat');

%  chanlocs = readlocs( 'd:\\Google_Drive\\MATLAB\\EEGLAB\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
%  load('D:\Google_Drive\MATLAB\LAB_MatlabScripts\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat');%ST-THETA-BURST
%  chanlocs=chanlocs62; %ST-THETA-BURST
for i=1: size(fileNames,1) %%%%%%%%%%%%%%%%% PIPELINE LOOP
%% loading the files
    
    fileName=fileNames{i};
    %fileName=Path(i).name;
    pathName=[Path(i).folder '\'];
% pathName='X:\TEMERTY\MDD_ECT_rTMS\rTMS\EEG_data\MDD266\PRE\';
% fileName='MDD266_PRE_TMS_DLPFC_SP2.cnt dec.cnt'
    if   strcmpi(fileName(:,end-2:end),'set')
        EEG = pop_loadset( [pathName fileName]);
    elseif strcmpi(fileName(:,end-2:end),'cnt')
% script for neuroscan CNT
        EEG = pop_loadcnt([pathName fileName] , 'dataformat', 'int32');
% script for ANT CNT
% EEG = pop_loadeep_v4([pathName fileName])
    elseif strcmpi(fileName(:,end-2:end),'cdt')
        EEG = loadcurry([pathName fileName], 'CurryLocations', 'True'); %'CurryLocations', 'False'
    elseif strcmpi(fileName(:,end-3:end),'vhdr')
        EEG  = pop_loadbv(pathName , fileName );
    end
% databak=EEG.data;
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

    EEG.session='1';%EEG.filename(19:22);
    EEG.data=double(EEG.data);
    EEG = eeg_checkset( EEG );
    % 
     disp(['Dataset is loaded ' num2str(i) '/' num2str(size(fileNames,1))])
% EEG = pop_saveset( EEG, [outdir EEG.setname]);
%% quality check plots
%  figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
%  pop_eegplot( EEG1, 1, 1, 1);
%  figure; pop_spectopo(EEG, 1, [-200  500], 'EEG' , 'percent', 50, 'freq', [6 10 22 45], 'freqrange',[2 60],'electrodes','off');
% figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'Oz'}),:)))
%% fix chanlocs - (change back to regular chanlocs!)
    %disp(['Fix Chanlocs  - dataset ' num2str(i) '/' num2str(size(fileNames,1))])
%   default TEMERTY cap
%   EEG=pop_chanedit(EEG, 'lookup','d:\\Google_Drive\\MATLAB\\EEGLAB\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');

%   EEG=pop_chanedit(EEG, 'lookup','A:\\WorkingSet\\\Pilot1\\CA-105.nlr_1-7-2022_10-20_AM.elc');   
%   figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
      
%   Shengze_tmseeg_202007_test_ZNN
%   EEG=pop_chanedit(EEG, 'lookup','A:\Workingset-A\shengze_tmseeg_202007_test_ZNN\ANT_64.ced');
     
       %chanBAK=EEG.chanlocs
    % sara trembley - theta burst bad motage fix
%     remove_channels={'TP9' 'TP10' 'M1' 'M2' 'VEO' 'HEO' 'EKG' 'EMG' 'HL 1' 'HL 2' 'Trigger'} ; % sara trembley - theta burst bad motage fix
%     EEG = NU_channel_removal_ST(EEG, remove_channels); % sara trembley - theta burst bad motage fix
%     chanlocs=EEG.chanlocs;
%  pop_eegplot( EEG, 1, 1, 1);    
    %EEG.chanlocs=chanlocs;
    %EEG=pop_chanedit(EEG, chanlocs); 
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo)
%% remove un-needed channels
    disp(['remove un-needed channels  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    remove_channels= {'lz' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo)
    % remove_channels= {'F11' 'F12' 'FT11' 'FT12' 'TP9' 'TP10' 'CB1' 'CB2' 'lz' 'M1' 'M2' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
    % remove_channels= {'TP9' 'TP10' 'CB1' 'CB2' 'lz' 'M1' 'M2' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
    chanlocs62=EEG.chanlocs;
        EEG = pop_select( EEG,'nochannel',remove_channels);
 
    % removing other channels that do not apear in the chanlocs files
   % EEG = pop_select( EEG,'nochannel',{EEG.chanlocs(~ismember({EEG.chanlocs.labels},{chanlocs62.labels})).labels});
%% Check for bridged electrodes (just for Quality Check)
%     disp(['check for bridged electrodes  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
%     EEG1=EEG;
%     
%     if EEG1.srate>1999
%         
%         rate=EEG1.srate/1000;
%         [~, ind_t0] = min(abs(EEG1.times-(0)));
%         T=EEG1.times; 
%         Dd=[smoothdata(EEG1.data(:,[ind_t0-1]:-1:1,:),2,'movmean',rate) smoothdata(EEG.data(:,ind_t0:end,:),2,'movmean',rate)]; %movmean gaussian
%         EEG1.times=[]; EEG1.times= [flip(T(ind_t0:-rate:1)) T(ind_t0+rate:rate:end)];
%         EEG1.data=[]; EEG1.data=[flip(Dd(:,ind_t0:-rate:1,:),2) Dd(:,ind_t0+rate:rate:end,:)];
%         EEG1.srate=EEG1.srate/rate;
%         EEG1.pnts=size(EEG1.data,2);
%         EEG1 = eeg_checkset( EEG1 );
%     clear Dd
%     end
%     [EB, ED]= eBridge(EEG1);
%     EEG.eBridge=EB.Bridged.Labels;
%     clear EEG1
%% Identify TMS events
    disp(['Identify TMS Pulse  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    %contains({EEG.event.type}
    evname= '128';   % evname= '1200001'
    newtrigs=1; %newtrigs=0
    method=1;
    DoublePulseINT=0; % paired-pulse interval, Zero (0) in case of single pulse
    [impotant, ~]=elecName(EEG,{'f5' 'f3' 'f1' 'fc5' 'fc3' 'fc1' 'fcz' 'fc2' 'c1' 'cz' 'c2' 'c4' 'c6' 'cp5' 'cp3'});
    if (contains(EEG.setname,'_RS'))
        [impotant, ~]=elecName(EEG,{'f6' 'f4' 'f2' 'fc6' 'fc4' 'fc1' 'fcz' 'fc2' 'c1' 'cz' 'c2' 'c4' 'c6' 'cp5' 'cp3'});
    end
% [impotant, ~]=elecName(EEG,{'f5' 'f3' 'fc5' 'fc3' 'fc4' 'fc6' 'f4' 'f6' 'cz' 'fcz'});
% [impotant ~]=elecName(EEG,{ 'fc5' 'fc3' });
% impotant=[7 8 9 16 17 18 ]; %region of stimulation 25 26 27
    EEG = IDpulse(EEG,method,DoublePulseINT,impotant,newtrigs,evname);
% pop_eegplot( EEG, 1, 1, 1);
% TESA PULSE DETECTION - TEST PULSE LABELED AS 'LICI', CONDITIONING PULSE LABELED AS 'con'  
% EEG.event=[];
% elec = 'cz';
% DoublePulseINT=101; % paired-pulse interval- LICI=101, single pulse = 0 
% 
% if DoublePulseINT==0
%     EEG = tesa_findpulsepeak(EEG, elec, 'dtrnd', 'poly', 'thrshtype','dynamic', 'wpeaks', 'gui', 'plots', 'off', 'tmsLabel', 'TMS','paired', 'no');
% else
%     EEG = tesa_findpulsepeak( EEG, elec, 'dtrnd', 'poly', 'thrshtype', 'dynamic', 'wpeaks', 'pos', 'plots', 'on', 'tmsLabel', 'TMS', 'paired', 'yes', 'ISI', DoublePulseINT, 'pairLabel', {'LICI'});
% end
% EEG = pop_saveset( EEG, 'filename',fileName,'filepath',['A:\WorkingSet\bipolar'],'check', 'on','savemode','onefile','version','7.3');
EEG = eeg_checkset( EEG );
%% Epoching (Wide)
    disp(['Epoching  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    % pop_eegplot( EEG, 1, 1, 1);
    %evname='Tms'; EEG1=EEG;
    ep=[-2.4  2.4]; % evname='128' evname='100128'  evname='230' ep=[-1.4  1.4]; ep=[-0.350 0.400];
    EEG = pop_epoch( EEG, {evname},ep , 'newname',EEG.setname );
    %   EEG = pop_epoch( EEG, {'R128'}, [-1.3 1.3]);
    EEG = Z_append(EEG,'_Epo');
    EEG = eeg_checkset( EEG );
    if size(size(EEG.data),2)<3
        error('epoching unsucsseful')
    return
    else
        disp(['Epoching  - sucsessful ' num2str(i) '/' num2str(size(fileNames,1))]);
    end
    mkdir([outdir '\EPO'])
    %EEG = pop_saveset( EEG, [ '111.set']);
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'EPO\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');

%%
    % figure; plot(EEG.times,squeeze(mean(EEG.data(elecName(EEG,{'f3'}),:,[50:80]),3)))
     % figure; plot(EEG.times,squeeze(mean(EEG.data(elecName(EEG,{'f3'}),:,1:100),3)))
     %figure; plot(EEG.times,squeeze(mean(EEG.data(elecName(EEG,{'cz'}),:,1:194),3)))
     % [~, c1] = min(abs(EEG.times-(-50)));[~, c2] = min(abs(EEG.times-(301)));
%
% figure; plot(EEG.times(c1:c2),squeeze(mean(EEG.data(elecName(EEG,{'f3'}),c1:c2,1:20),3)))
     % figure; plot(EEG.times([c1:c2]),squeeze(EEG.data(elecName(EEG,{'af3'}),c1:c2,1:100)))
       % pop_eegplot( EEG, 1, 1, 1);
       % figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'c4'}),:,:)))
         %figure; plot(EEG.times,squeeze(mean(EEG.data(elecName(EEG,{'f3'}),:,:),3)))
          %figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:)))
          %%%
%           [x,m,n]=size(EEG.data);
%accumarray([1:5:n],EEG.data(elecName(EEG,{'f3'}),:,:)))
%           dat=reshape(EEG.data(elecName(EEG,{'f3'}),:,:),1,[m*n])
%           edges = (0:5:length(dat));
% [~,~,loc]=histcounts(EEG.times,edges);
% meany = accumarray(loc(:),y(:))./accumarray(loc(:),1);
% xmid = 0.5*(edges(1:end-1)+edges(2:end));
%     mkdir([outdir '\Epoched'])
%     %EEG = pop_saveset( EEG, [outdir '\Epoched\' EEG.filename]);
%     EEG = pop_saveset( EEG, 'filename','222.set','check', 'on','savemode','onefile','version','7.3');

%     EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'Epoched\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
%% Trim the TMS-pulse 
    disp(['Trim the TMS-pulse  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    trim=[-2 15]; %DoublePulseINT=0 % trim=[-2 5];
    EEG = Triminterp(EEG,trim);
    if DoublePulseINT>0
        EEG = Triminterp(EEG,trim-DoublePulseINT);
    end
   
    %pop_eegplot( EEG, 1, 1, 1);
    %mkdir([outdir '\Trimed'])
    %EEG = pop_saveset( EEG, [outdir '\DEC\' EEG.filename]);
    %EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'Trimed\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
%% DeTrending Epochs (also baseline correct)
    disp(['Detrending Epochs  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    % EEG = pop_loadset( ['D:\WORDKINGset_D\Optimal_targetting_SP\Epoched\OT_017_OJB3_1MV_SP_MASKING_Epo.set']);
    % EEG1.data=EEG.data;
    EEG.data=double(EEG.data);
    %disreg=[-DoublePulseINT-350 350];
    %disreg_vec=EEG.times<disreg(1) | EEG.times>disreg(2); 
[~, bc1] = min(abs(EEG.times-(-160)));
[~, bc2] = min(abs(EEG.times-(-20)));
%     [~, bp2] = min(abs(EEG.times-(disreg(2))));
%     %disreg_vec=[[1:20:bp1] [bp2:20:size(EEG.times,2)]] ;
%     [~, bp3] = min(abs(EEG.times));
%     breakpoint=[bp1 bp3 bp2 ];
    for elecc=1:size(EEG.data,1)
        for epoo=1:size(EEG.data,3)
            EEG.data(elecc,:,epoo)=detrend(double(EEG.data(elecc,:,epoo)),1);
            EEG.data(elecc,:,epoo)=EEG.data(elecc,:,epoo)-mean(EEG.data(elecc,[bc1:bc2],epoo),'omitnan'); % baseline correction
%              figure;plot(EEG.times,double(EEG.data(elecc,:,epoo)));
%              yoyoyo= EEG.data(elecc,:,epoo);
%              plot(EEG.times,detrend(double(EEG.data(elecc,:,epoo)),1,breakpoint,'Continuous',true),'-r')
%                plot(EEG.times,detrend(double(EEG.data(elecc,:,epoo)),1,'Continuous',true),'-r')
%              figure;plot(EEG.times(1:breakpoint),double(EEG.data(elecc,1:breakpoint,epoo))); 
%              plot(EEG.times(1:breakpoint),detrend(double(EEG.data(elecc,1:breakpoint,epoo)),1),'-r')
            %EEG.data(elecc,:,epoo)=detrend(double(EEG.data(elecc,:,epoo)),2,'Continuous',0,'SamplePoints',disreg_vec);
            %EEG.data(elecc,:,epoo)=detrend(double(EEG.data(elecc,:,epoo)),2,disreg_vec);
        end
    end
    EEG = eeg_checkset( EEG );        
    %figure; hold on,  plot(EEG.times,EEG.data(elecc,:,epoo))
    %plot(EEG.times,detrend(double(EEG.data(elecc,:,epoo)),2,disreg_vec),'-r')
    %plot(detrend(EEG.data(elecc,:,epoo),1))
    %plot(EEG.times,EEG.data(elecc,:,epoo))
    %pop_eegplot( EEG, 1, 1, 1);
   % mkdir([outdir '\detrend'])
    %EEG = pop_saveset( EEG, [outdir '\DEC\' EEG.filename]);
   % EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'detrend\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
%% remove bad channels
disp(['remove bad channels  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);

imp={'fp1' 'fp2' 'fpz' 'fz' 'fc3' 'fc5' 'f3' 'f5'};
if (contains(EEG.setname,'_RS'))
    imp={'fp1' 'fp2' 'fpz' 'fz' 'fc4' 'fc6' 'f4' 'f6'};
end
% imp={'fp1' 'fp2' 'fpz' 'c5' 'c3'}; this line is for motor stimulation
% see if you can introduce trial by trial eeglab plugin
if EEG.srate>1999
EEG1=resample_NOfilt(EEG, 1000);
end

if contains(lower(fileName),'_pp')
    disreg=[-DoublePulseINT-5 400];
elseif contains(lower(fileName),'_sp')
    disreg=[-5 400];
else
    disreg=[-5 400];
end

[~, t1] = min(abs(EEG1.times-(disreg(1))));
[~, t2] = min(abs(EEG1.times-(disreg(2))));
%EEG1.times(t1:t2)=zeros(1,size(EEG1.times(t1:t2),2));
EEG1.data(:,t1:t2,:)=zeros(size(EEG1.data,1),size(EEG1.data(:,t1:t2,:),2),size(EEG1.data,3));
EEG1.data=double(EEG1.data);
%EEG1.pnts=length(EEG1.times);
EEG1 = eeg_checkset( EEG1 );
%pop_eegplot( EEG1, 1, 1, 1);
%withOUT important channels

EEG1.rejelecs={};
% remove bad channels: dead chennels
% check std of entire recording DROP=[2 3 4]
[DROP,  dropNaN, dropdead, drophigh] = dropOFFchan( EEG1 );
%EEG1.rejelecs=[EEG1.chanlocs([DROP]).labels];

if  isfield(EEG.chanlocs,'median_impedance')
    EEG1.rejelecs = unique( {EEG1.chanlocs([DROP]).labels EEG1.chanlocs(logical([EEG.chanlocs.median_impedance]>25)).labels} ) ;
    EEG1.rejelecs_impedence ={EEG1.chanlocs(logical([EEG.chanlocs.median_impedance]>25)).labels}
else
    EEG1.rejelecs = unique( {EEG1.chanlocs([DROP]).labels} ) ;
end
%     impchan=ismember(lower(imp),lower(EEG1.rejelecs));
%     EEG.rejimpchan=imp(impchan);
%     EEG.rejElecs=EEG1.rejelecs;
%EEG1 = pop_select( EEG1,'nochannel',EEG1.rejelecs);

important=elecName(EEG1,imp);
noimp=ones(1,size({EEG1.chanlocs.labels},2)); noimp(1,important)=0; noimp=logical(noimp); a=1:size(EEG1.data,1);aa= a(noimp);
[~, rejelecs] = pop_rejchan(EEG1, 'elec',aa ,'threshold',2.5,'norm','on','measure','spec','freqrange',[0.3 6] );
EEG1.rejelecs = unique( [{EEG1.chanlocs([aa(rejelecs)]).labels} EEG1.rejelecs] ) ;

important=elecName(EEG1,imp);
noimp=ones(1,size({EEG1.chanlocs.labels},2)); noimp(1,important)=0; noimp=logical(noimp); a=1:size(EEG1.data,1);aa= a(noimp);
[~, rejelecs] = pop_rejchan(EEG1, 'elec',aa ,'threshold',2.5,'norm','on','measure','spec','freqrange',[58 62] );
EEG1.rejelecs = unique( [{EEG1.chanlocs([aa(rejelecs)]).labels} EEG1.rejelecs] ) ;

%with important channels
[~, rejelecs] = pop_rejchan(EEG1, 'elec',[1:size(EEG1.data,1)] ,'threshold',3.5,'norm','on','measure','spec','freqrange',[0.3 6] );
EEG1.rejelecs = unique( [{EEG1.chanlocs([rejelecs]).labels} EEG1.rejelecs] ) ;

[~, rejelecs] = pop_rejchan(EEG1, 'elec',[1:size(EEG1.data,1)] ,'threshold',3.5,'norm','on','measure','spec','freqrange',[58 62] );
EEG1.rejelecs = unique( [{EEG1.chanlocs([rejelecs]).labels} EEG1.rejelecs] ) ;

[~, rejelecs] = pop_rejchan(EEG1, 'elec',[1:size(EEG1.data,1)] ,'threshold',6,'norm','on','measure','prob');
EEG1.rejelecs = unique( [{EEG1.chanlocs([rejelecs]).labels} EEG1.rejelecs] ) ;


%     impchan=ismember(lower(imp),lower(EEG1.rejelecs));
%     EEG.rejimpchan=imp(impchan);
%     EEG.rejElecs=EEG1.rejelecs;
EEG = pop_select( EEG,'nochannel',EEG1.rejelecs);
EEG.rejelecs=EEG1.rejelecs;
%pop_eegplot( EEG, 1, 1, 1);
[EEG, rejepochs] = pop_autorej(EEG, 'nogui','on','eegplot','off');
EEG.rejEpochs=rejepochs;

clear EEG1
%% downsampling
    disp(['downsampling  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    
    if EEG.srate>1999
    %EEG = pop_resample(EEG, 1000);     
    EEG=resample_NOfilt(EEG, 1000);
    end

    %    pop_eegplot( EEG, 1, 1, 1, {}, 'spacing', 70)
    mkdir([outdir '\DEC'])
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'DEC\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
%% Channel Interpolation - Better not to run at this stage
%     disp(['Channel Interpolation  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
%     EEG = pop_interp(EEG, chanlocs62, 'spherical');
% %   EEG = pop_chanedit(EEG, 'lookup','E:\\Google_Drive\\MATLAB\\EEGLAB\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc');
%     remove_channels2= {'VEO' 'HEO' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'}; %'TP9' 'TP10' 'CB1' 'CB2' 'M1' 'M2' 
%     EEG = pop_select( EEG,'nochannel',remove_channels2);
%     EEG = Z_append(EEG,'_interp');
%     EEG = eeg_checkset( EEG );
%     mkdir([outdir '\interp'])
%     EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'interp\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
%     % EEG = pop_saveset( EEG, 'filename','111.set','filepath',strrep([outdir 'interp\'],'\','\\'));
%     %EEG = pop_saveset( EEG, 'filename','111.set','filepath','D:\\WORDKINGset_D\\Optimal_targetting_SP\\interp\\');
%     % EEG = pop_loadset( ['D:\WORDKINGset_D\Optimal_targetting_SP\interp\OBJ3_OT18_1mv_SP_MASKING_Epo_Dec_interp.set']);
%     %pop_eegplot( EEG, 1, 1, 1);
%% Baseline correction
    disp(['Baseline correction  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    base=[-200   -20]; %  base=[-200   -20];
    EEG  = pop_rmbase( EEG, base );
    if DoublePulseINT>0
        EEG  = pop_rmbase( EEG, base-DoublePulseINT);
    end
    EEG  = eeg_checkset( EEG );
    %     
    % EEG_test=EEG;
    % %if ~strcmp(EEG_test.ref,'averef'), EEG_test = pop_reref( EEG_test, []); end
    % ch_ROI=[8 ]; %ERP for region of intrest
    % figure,
    % plot(EEG_test.times, squeeze(mean(mean(EEG_test.data(ch_ROI,:,:),3),1)))
%% ICA1
    disp(['ICA 1  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    EEG.data=double(EEG.data);
    %EEG = pop_runica(EEG, 'lrate', 0.001, 'extended',1, 'interupt','off','verbose','off');
    %EEG1=EEG;
    % check wICA
    %[wIC, A, W, IC] = wICA(EEG,'runica', 1, 0, srate, 5); slower but better 
 %    [wIC, A, W, IC] = wICA(EEG,'fastica', 1, 0, srate, 5);
 %     artifacts = A*wIC;
     %reshape EEG signal from EEGlab format to channelsxsamples format
 %   EEG2D=reshape(EEG.data, size(EEG.data,1), []);
    %subtract out wavelet artifact signal from EEG signal
 %   wavcleanEEG=EEG2D-artifacts;
 %   EEG.data = wavcleanEEG;
 %EEG = pop_runica(EEG, 'icatype','fastica','verbose','off');
    %EEG = pop_runica(EEG,'icatype','binica', 'extended',1,'interupt','on','pca',size(EEG.data,1) ); % 
    %EEG = pop_runica(EEG,'icatype','fastica','approach','symm','g','tanh','stabilization','on');
    EEG = pop_runica(EEG, 'icatype', 'picard','m',12 , 'maxiter',700); %,'mode','standard';
    %save data
    % create ICA folder
    EEG = Z_append(EEG,'_ICA1');
    mkdir([outdir '\ICA1'])
    %EEG = pop_saveset( EEG, [outdir 'ICA1\' EEG.filename]);
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'ICA1\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
%% removes 'decay' & line-noise ICA components
    disp(['ICA1 - Comps removal  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    
    
    if (contains(EEG.setname,'_RSP'))
        ROI = {'F4' 'fpz' 'fp2' 'AF4' 'F12' 'FT12' 'F2' 'FT8' 'AF6' 'F8' 'FC8' 'FC6' 'F6' 'FC4' 'c4' 'c6' 't8' 'af8' 'afz' 'fz' 'fcz' 'FP2'};
    elseif (contains(EEG.setname,'_LSP'))
        ROI = {'F3' 'fpz' 'fp2' 'AF3' 'F11' 'FT11' 'F1' 'FT7' 'AF5' 'F7' 'FC7' 'FC5' 'F5' 'FC3' 'c3' 'c5' 't7' 'af7' 'afz' 'fz' 'fcz' 'FP1'};
    end

    EEG = ICAremDECAY(EEG,ROI);
    minfreq=58; maxfreq=62;
    EEG = FFT_comp(EEG, minfreq, maxfreq);
 
   % figure; plot(EEG.times,squeeze(mean(EEG.data(elecName(EEG,{'f3'}),:,[50:80]),3)))
     % figure; plot(EEG.times,squeeze(mean(EEG.data(elecName(EEG,{'f3'}),:,1:100),3)))
     % [~, c1] = min(abs(EEG.times-(-50)));[~, c2] = min(abs(EEG.times-(301)));%
% figure; plot(EEG.times(c1:c2),squeeze(mean(EEG.data(elecName(EEG,{'f3'}),c1:c2,:),3)))
%% Filtering after removing the decay artifact
    disp(['Filtering  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    EEG.data=double(EEG.data);
    %  pop_eegplot( EEG, 1, 1, 1);
    low=0.2;
    high=63;
    %ACTIVE_New_2WayFilter(1,48,0)
    EEG.data = FilterTEP(EEG.data, low, high,EEG.srate); % ,[1:size(EEG.data,1)],EEG.srate
    % FilterTEP(EEGdata, low, high,elecs,srate)
     %     pop_comperp( EEG, 1, [1],[],'tlim',[1:450],'chans',[ 1:size(EEG.chanlocs,2)],'addavg','off','addstd','off','addall','on','diffavg','off','diffstd','off','tplotopt',{'ydir' 1});

%plot(EEG.times,EEG.data(15,:,20))
    %EEG = Z_append(EEG,'_pruned1_filt');
    %EEG = pop_saveset( EEG, [ EEG.filename]);
    % save data *_rmDecay_Filterd1_48.set
%% re-referencing
    % disp(['channel re-referencing  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    % EEG.data=double(EEG.data);
    % EEG = pop_reref(EEG, []);
%% DeTrending Epochs Second time
    disp(['Detrending Epochs  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    EEG.data=double(EEG.data);
    for elecc=1:size(EEG.data,1)
        for epoo=1:size(EEG.data,3)
            EEG.data(elecc,:,epoo)=detrend(double(EEG.data(elecc,:,epoo)),1);
        end
    end
    EEG = eeg_checkset( EEG );        
    
    %pop_eegplot( EEG, 1, 1, 1);
%% Baseline correction
    disp(['Baseline correction  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    base=[-800   -200];
    EEG  = pop_rmbase( EEG, base );
    if DoublePulseINT>0
        EEG  = pop_rmbase( EEG, base-DoublePulseINT);
    end
    EEG  = eeg_checkset( EEG );
    
    %    pop_eegplot( EEG, 1, 1, 1, {}, 'spacing', 70)
%% ICA2
    disp(['ICA2  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    EEG.data=double(EEG.data);
    %EEG = pop_runica(EEG, 'lrate', 0.0001, 'extended',1, 'interupt','off','verbose','off');
    %EEG = pop_runica(EEG,'icatype','fastica','approach','symm','g','tanh','stabilization','on');.
    EEG = pop_runica(EEG, 'icatype', 'picard','m',12 , 'maxiter',700); %,'mode','standard';
   % EEG = pop_runica(EEG,'icatype','binica','extended',1,'verbose','off')
    EEG = Z_append(EEG,'_ICA2');
    mkdir([outdir '\ICA2'])
   % EEG = pop_saveset( EEG, [ 'ICA2\' EEG.filename]);
   EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([outdir 'ICA2\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
     %    pop_eegplot( EEG, 1, 1, 1, {}, 'spacing', 70)
%% Clear
    clearvars -except chanlocs i fileNames pathName Path outdir
end 
%% FIGURES - QUALITY CHECK
% pathName='D:\Google_Drive\PhD Zangen\ADHD\EEGLAB GRANDS BACKUP\TEP - 15ms\';
% fileNames=[{'ADHD GRAND Merged datasets Unstich-basecorr.set'}, ...
%     {'Healthy GRAND Merged datasets Unstich-basecorr.set'}];
% ALLEEG = pop_loadset('filepath' , pathName, 'filename',fileNames);
 % pop_eegplot( EEG, 1, 1, 1);

time= [-30 400]
elec='F3'
for dset=1:8
[~, T1] = min(abs(ALLEEG(dset).times-(time(1))));
[~, T2] = min(abs(ALLEEG(dset).times-(time(2))));
% plot channel trialXtimeXmvolt heatmap
%figure; pop_erpimage(EEG,1, elecName(EEG,{elec}),[[]],elec,10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [8] EEG.chanlocs EEG.chaninfo } );
% plot one channel traces
% figure; hold on; plot(ALLEEG(dset).times(T1:T2),squeeze(ALLEEG(dset).data(elecName(ALLEEG(dset),{elec}),T1:T2,:))); 
% yy1=min(squeeze(ALLEEG(dset).data(elecName(ALLEEG(dset),{elec}),T1:T2,:)),[],'all')-5;
% yy2=max(squeeze(ALLEEG(dset).data(elecName(ALLEEG(dset),{elec}),T1:T2,:)),[],'all')+5;
% fill([0 15 15 0],[yy2 yy2 yy1 yy1],'w','edgecolor','none');  
% title([strrep(ALLEEG(dset).setname,'_',' ') ' ' elec])
% hold off
% plot actication headplot
figure; pop_timtopo(ALLEEG(dset), time, [NaN], ALLEEG(dset).setname,'title',[strrep(ALLEEG(dset).setname,'_',' ')] ); 
%   figfig('F3')

end
ALLEEG=EEG;
figure; pop_timtopo(ALLEEG(dset), time, [NaN] )
time= [-30 150]
pop_comperp( ALLEEG, 1, [4 5 6 7 8 9 10 11] ,[],'addavg','off','addstd','off','addall','on','diffavg','off','diffstd','off','tplotopt',{'ydir' 1},'tlim',time );
pop_comperp( ALLEEG, 1, [1 2 3 4 5 6 7 8] ,[],'addavg','off','addstd','off','addall','on','diffavg','off','diffstd','off','tplotopt',{'ydir' 1},'tlim',time );
pop_comperp( ALLEEG, 1, [1 2 3] ,[],'addavg','off','addstd','off','addall','on','diffavg','off','diffstd','off','tplotopt',{'ydir' 1},'tlim',time );

%figfig(
%% POST-processing
%% Concatenate the experimental conditions datasets
%     load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat');
%     EEG = pop_interp(EEG, chanlocs62, 'spherical');
%    subs=
%    cond=
%     errrorsss = concatenate(pathName, subs, cond);
%           
%% load
clear all
pathi='D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\ICA2\';
[fileNames, pathName]=Z_getSetsFileNames('set',pathi);
i=0
%%
i=i+1; %hhh
fileName=fileNames{i};
EEG = pop_loadset( [pathName fileName]);
[EEG, ~] = tesa_sortcomps(EEG);
EEG = pop_iclabel(EEG, 'default');
RejComp=sort([find(any(EEG.etc.ic_classification.ICLabel.classifications(:,4:6)>0.9,2));...
    find(any(EEG.etc.ic_classification.ICLabel.classifications(:,2:3)>0.86,2))]);
%EEG = pop_subcomp( EEG,RejComp , 0);

%% TEP before and after component removal
%[EEG, ~] = tesa_sortcomps(EEG);
%TMS_pop_selectcomps(EEG, [1:35] );
%TMS_pop_selectcomps(EEG, [36:55] );
% pop_selectcomps(EEG, [1:35] );
% pop_selectcomps(EEG, [36:55] );
EEG = pop_iclabel(EEG, 'default');
EEG = pop_viewprops( EEG, 0, [1:35], {'freqrange', [2 60]}, {}, 1, '' )
% EEG = pop_select( EEG,'nochannel',7);
% pop_eegplot( EEG, 1, 1, 1); 
% plot(EEG.times,EEG.data(5,:,2))
pop_plotdata(EEG, 0, [1:size(EEG.icaweights,1)] , [1:size(EEG.data,3)], EEG.setname, 0, 1, [-3 8]);
%% Before-After ICA copmponent removal
%[EEG, varsPerc] = tesa_sortcomps(EEG);
channels= elecName(EEG,{'f3'}) %41  27 47 20
%components= [ 1 2 ]
%components= [find(EEG.reject.gcompreject)] %3     4     5     6     7     8     9    11    35
components= RejComp'
comp='on';
comprem(EEG,channels,components,comp)
% pop_eegplot( EEG, 1, 1, 1, {}, 'spacing', 70)
%iclabel
%% save file after compnents removal
%components= [find(EEG.reject.gcompreject)];
pathh='D:\WORKINGset-D\ECT\PRE_SP\PRE_ECT_cleanDATA\'; % CHANGE PATH!!!!

if isfield(EEG,'Removed_ICA_Comps')
if isfield(EEG.Removed_ICA_Comps,'remcomp_ICA2')
    EEG.Removed_ICA_Comps.remcomp_ICA2_1=components;
else
    EEG.Removed_ICA_Comps.remcomp_ICA2=components;
end
else
    EEG.Removed_ICA_Comps.remcomp_ICA2=components;
end

% ---- Channel Interpolation 
     disp(['Channel Interpolation  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
     EEG = pop_interp(EEG, chanlocs, 'spherical');
% %   EEG = pop_chanedit(EEG, 'lookup','E:\\Google_Drive\\MATLAB\\EEGLAB\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc');
     remove_channels2= {'CB1' 'CB2' 'VEO' 'HEO' 'EKG' 'EMG' 'HL 1' 'HL 2' 'Trigger'};
     EEG = pop_select( EEG,'nochannel',remove_channels2);
     EEG = Z_append(EEG,'_interp');
     EEG = eeg_checkset( EEG );
%     mkdir([pathName '\interp'])
%     EEG = pop_saveset( EEG, [pathName '\interp\' EEG.filename]);
%     %pop_eegplot( EEG, 1, 1, 1);

setname1=EEG.setname;
filename1=EEG.filename;
datfile1=EEG.datfile;
EEG = pop_subcomp( EEG, components, 0);
EEG.setname=setname1;
EEG.filename=filename1;
EEG.datfile=datfile1;
EEG=Z_append(EEG, '_CLEEG');
EEG = eeg_checkset( EEG );

EEG = pop_saveset( EEG, [pathh EEG.filename]);
if CURRENTSET==size(ALLEEG,2)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',1,'study',0);   
else
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',CURRENTSET+1,'study',0); 
end

[EEG, varsPerc] = tesa_sortcomps(EEG);
TMS_pop_selectcomps(EEG, [1:35] );
pop_plotdata(EEG, 0, [1:size(EEG.icaweights,1)] , [1:size(EEG.data,3)], EEG.setname, 0, 1, [-3 8]);
disp('Current dataset'), disp(CURRENTSET);
%% GRAND Average 
%chan_interp='off';
%load('D:\Google_Drive\MATLAB\LAB_MatlabScripts\chanlocs68.mat');
%load('D:\MATLAB\LAB_MatlabScripts\Chanlocs\chanlocs66_flexnet_compumedics.mat'); chanlocs=chanlocs62;
%load('D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\channel_locations\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat'); chanlocs=chanlocs62;
pathi='D:\DATA\WellcomeLeap_TMS-EEG\RAW_SP\interp\ICA2\final1\'
chan_interp='on'
chanlocs=chanlocs66;
STDcalc=0
Z2_grand_average('WL_SP_all',[1:6],pathi, chan_interp,chanlocs,STDcalc) 

figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:))); ylim([-40 40]); xlim([-300 300])
figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:))); ylim([-40 40]); xlim([-300 300])


[a loc]=max(EEG.data(17,5,:))
EEG.data(:,:,8)=[]
%% SNR run over grandaverage

[~, bc1] = min(abs(EEG.times-(-180)));
[~, bc2] = min(abs(EEG.times-(-100)));
[~, bc3] = min(abs(EEG.times-(20)));
[~, bc4] = min(abs(EEG.times-(100)));
snr(squeeze(EEG.data(elecName(EEG,{'f3'}),[bc3:bc4],:)),squeeze(EEG.data(elecName(EEG,{'f3'}),[bc1:bc2],:)))
figure; hold on; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:))); ylim([-40 40]); xlim([-300 300])
fill([-180 -100 -100 -180],[-40 -40 40 40],'-k','facealpha',.15,'edgecolor','none'); 
fill([20 100 100 20],[-40 -40 40 40],'-k','facealpha',.15,'edgecolor','none'); 
fill([-2 20 20 -2],[-40 -40 40 40],'w','edgecolor','none'); 

figfig('signal2noise')
%% old _Auto ICA-comp removal
Z_TESA_ICAautoeyes

BADchan_comp
%%

%plotting grand average per subject time-course
figure; hold on;
for j=1:size(EEG.data,3)
plot(EEG.times,EEG.data(28,:,j))
end

timel=[140 215];
[~, t1] = min(abs(EEG.times-(timel(1))));
[~, t2] = min(abs(EEG.times-(timel(2))));
 [mm subs]=max(mean(EEG.data(28,t1:t2,:),2))
 
 figure('position', [100 300 1500 500]); hold on;
  for jj=1:size(EEG.data,1)
 plot(EEG.times(1700:2400),mean(EEG.data(jj,1700:2400,:),3))
 end
 title(ALLEEG(CURRENTSET).filename)
 export_fig LICI.pdf -q101 -painters -append
%% paired-Pulse GRANDaverage
Z2_grand_average('non_jitter_1_ICA',[4:6]) 

newsetname='LICI'
int=100
Z_Grand_diff_paired(newsetname,int)

figure;plot(EEG.times,squeeze(EEG.data(17,:,:)))
%% deconcatenate concatenate
pathName='D:\WORKINGset-D\ECT\PRE_SP\removed_comps_ICA2\concatenate\';
subloc=[4:6];
condloc=[22:24];
concatenate(pathName, subloc , condloc)


deconcatenate
%% WAVELET !
%TMS_NewWavelet2
TMS_NewWavelet3

%%
EEG=Z_append(EEG, '_CLEEG');



EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, [pathName fileName]);
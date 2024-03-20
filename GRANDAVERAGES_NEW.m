function GRANDAVERAGES_NEW(EventName, Subsind,pathi,fileNames,chan_interp,chanlocs)
    eeglabdir='/media/ipp/DATA/Documents/MATLAB/eeglab2024.0';

for i=1: size(fileNames,1)
    fileName=fileNames{i};
    EEG = pop_loadset( [pathi fileName]);
    if strcmp(chan_interp,'on')
        EEG = pop_interp(EEG, chanlocs, 'spherical');
        remove_channels= { 'F11' 'F12' 'FT11' 'FT12' 'CB1' 'CB2'...
        'lz' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
        EEG = pop_select( EEG,'nochannel',remove_channels);
    end
    
    if i==1
        GrandAverage = zeros(size(EEG.data, 1), size(EEG.data, 2), size(fileNames,1));
               
        if isfield(EEG,'preprosTAB')
         preprosTAB=table;
         preprosTAB=EEG.preprosTAB;
        end
    end
    
    
    GrandAverage(:,:,i) =nanmean(EEG.data, 3); %  nanmean(EEG.data, 3);
        
    if isfield(EEG,'preprosTAB') & i>1
    preprosTAB=[preprosTAB; EEG.preprosTAB];
    else
    subjects(1, i)={[[fileName(Subsind)]]}; %{['3' [fileName(Subsind)]]}
    end
    %%
    METADATA{:,i}.Subject=EEG.subject;
    METADATA{:,i}.filename=EEG.setname;
    METADATA{:,i}.group=EEG.group;
    METADATA{:,i}.condition=EEG.condition;
    METADATA{:,i}.session=EEG.session;
    METADATA{:,i}.trials=EEG.trials;
    METADATA{:,i}.srate=EEG.srate;
    METADATA{:,i}.epoch_num=size(EEG.data,3);
    METADATA{:,i}.event=EEG.event;
    METADATA{:,i}.urevent=EEG.urevent;
    METADATA{:,i}.epochs=EEG.epoch;
    if isfield(EEG,'eBridge')
        METADATA{:,i}.eBridge=EEG.eBridge;
    end
    if isfield(EEG,'tmscut')
        METADATA{:,i}.tmscut=EEG.tmscut;
    end
    if isfield(EEG,'rejelecs')
        METADATA{:,i}.rejelecs=EEG.rejelecs;
    end
    if isfield(EEG,'rejEpochs')
        METADATA{:,i}.rejEpochs=EEG.rejEpochs;
    end
    if isfield(EEG,'Removed_ICA_Comps')
        if  isfield(EEG.Removed_ICA_Comps,'comps_num')
            METADATA{:,i}.comps_num=EEG.Removed_ICA_Comps.comps_num;
        end
        if  isfield(EEG.Removed_ICA_Comps,'remcomp_ICA2')
            METADATA{:,i}.remcomp_ICA2=EEG.Removed_ICA_Comps.remcomp_ICA2;
        end
    end
    %%
end


EEG.data = GrandAverage;
EEG.subject='GRAND';
EEG.condition=EventName; %(max(Subsind)+2:max(Subsind)+5);

if isfield(EEG,'preprosTAB') & i>1
    EEG.preprosTAB=table;
    EEG.preprosTAB=preprosTAB;
else
    EEG.subjects=subjects;
end

EEG.METADATA=METADATA;
EEG.icaact=[];
EEG.icawinv=[];
EEG.icachansind=[];
EEG.icasphere=[];
EEG.icaweights=[];
EEG.icachansind=[];
EEG.icasplinefile=[];
EEG.Removed_ICA_Comps=[];
EEG.rejEpochs=[];
EEG.rejelecs=[];
EEG.eBridge=[];


EEG.trials = size(fileNames,1);

EEG.event(find(strcmpi({EEG.event.type},'old')))=[];
for ee=1:size(fileNames,1)
    EEG.event(ee).type='Subject';
    EEG.event(ee).latency=[size(EEG.times,2).*ee+1];
    EEG.event(ee).duration=1;
    EEG.urevent(ee).channel=[];
    EEG.urevent(ee).bvtime=[];
    EEG.urevent(ee).bvmknum=ee;
    EEG.urevent(ee).code='Response';
    EEG.event(ee).urevent=ee;
    EEG.event(ee).epoch=ee;
    EEG.urevent(ee).latency=size(EEG.times,2).*ee+1;
    EEG.urevent(ee).duration=1;
    EEG.urevent(ee).channel=0;
    EEG.urevent(ee).bvtime=[];
    EEG.urevent(ee).bvmknum=ee;
    EEG.urevent(ee).type='Subject';
    EEG.urevent(ee).code='Response';
    EEG.epoch(ee).event=ee;
    EEG.epoch(ee).eventlatency=0;
    EEG.epoch(ee).eventduration=0;
    EEG.epoch(ee).eventchannel=0;
    EEG.epoch(ee).eventbvtime=[];
    EEG.epoch(ee).eventbvmknum=ee;
    EEG.epoch(ee).eventtype='Subject';
    EEG.epoch(ee).eventcode='Response';
    EEG.epoch(ee).eventurevent=ee;
end

EEG.event=EEG.event(1:size(fileNames,1));
EEG.urevent=EEG.urevent(1:size(fileNames,1));
EEG.epoch=EEG.epoch(1:size(fileNames,1));


EEG.setname=[EventName '_GrandAverage'];
EEG.filename=[EEG.setname '.set'];
EEG.datfile=[EEG.setname '.fdt'];
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, [pathi EEG.filename]);


return;


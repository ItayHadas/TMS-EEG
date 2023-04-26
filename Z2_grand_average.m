function Z2_grand_average(EventName, Subsind,chan_interp,chanlocs,STDcalc)
%EventName='WL_SP_PRE'; Subsind=[1:6],chan_interp,chanlocs,STDcalc=0
% EventName='ThetaB_1_Pre',Subsind=[2:10]
% EventName='ADHD1_1-15ms-ala_rogash'
% EventName='MST_motor_HC'   EventName='pre_MST2a'
% EventName='post_rtms'
% Subsind=4:6
[fileNames, pathName]=Z_getSetsFileNames('set');
% EventName='Healthy1_1-exponent_15ms'; Subsind=[5:7];
% EventName='LICI1_non'
for i=1: size(fileNames,1)
    
    if size(fileNames, 1)==1
        fileName=fileNames{i,1}';
    else
        fileName=fileNames{i,1};
    end
    
    EEG = pop_loadset( [pathName fileName]);
    [ num2str(i) ' out of ' num2str(size(fileNames,1)) '   -   ' fileName ]
    
    if strcmp(chan_interp,'on')
        EEG = pop_interp(EEG, chanlocs, 'spherical');
        %   EEG = pop_chanedit(EEG, 'lookup','E:\\Google_Drive\\MATLAB\\EEGLAB\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc');
        remove_channels2= {'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'HL 1' 'HL 2' 'Trigger'};
        EEG = pop_select( EEG,'nochannel',remove_channels2);
    end
    
    if i==1
        %GrandAverage = zeros(size(EEG.ERP, 1), size(EEG.ERP, 2), size(fileNames,1));
        %GrandAverage = nan(size(EEG.ERP, 1), size(EEG.ERP, 2), size(fileNames,1));
        GrandAverage = zeros(size(EEG.data, 1), size(EEG.data, 2), size(fileNames,1));
        if STDcalc==1
        GrandSTD = zeros(size(EEG.data, 1), size(EEG.data, 2), size(fileNames,1));
        end
        
        if isfield(EEG,'preprosTAB')
         preprosTAB=table;
         preprosTAB=EEG.preprosTAB;
        end
    end
    
    
    GrandAverage(:,:,i) =nanmean(EEG.data, 3); %  nanmean(EEG.data, 3);
    if STDcalc==1
    GrandSTD(:,:,i) =nanstd(EEG.data,0, 3); %
    end
    
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
if STDcalc==1
EEG.STDdata = GrandSTD;
end
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
%EEG = pop_newset( EEG, 1,'savenew',[EEG.filename],'gui','off');

EEG = pop_saveset( EEG, [pathName EEG.filename]);


return;

%%
 
 % EEG.METADATA=table;
% for s=1:size(METADATA,2)
% EEG.METADATA(s,:)=struct2table(METADATA{s});
% end



    %         METADATA{:,i}.Subject=EEG.subject;
    %         METADATA{:,i}.filename=EEG.setname;
    %         METADATA{:,i}.group=EEG.group;
    %         METADATA{:,i}.condition=EEG.condition;
    %         METADATA{:,i}.session=EEG.session;
    %         METADATA{:,ii}.trials=EEG.trials;
    %         METADATA{:,i}.srate=EEG.srate;
    %         METADATA{:,i}.epoch_num=size(EEG.data,3);
    %         METADATA{:,i}.event=EEG.event;
    %         METADATA{:,i}.urevent=EEG.urevent;
    %         METADATA{:,i}.epochs=EEG.epoch;
    %         if isfield(EEG,'eBridge')
    %         METADATA{:,i}.eBridge=EEG.eBridge;
    %         end
    %         if isfield(EEG,'tmscut')
    %         METADATA{:,i}.tmscut=EEG.tmscut;
    %         end
    %         if isfield(EEG,'rejelecs')
    %         METADATA{:,i}.rejelecs=EEG.rejelecs;
    %         end
    %         if isfield(EEG,'rejEpochs')
    %         METADATA{:,i}.rejEpochs=EEG.rejEpochs;
    %         end
    %         if isfield(EEG,'Removed_ICA_Comps')
    %           if  isfield(EEG.Removed_ICA_Comps,'comps_num')
    %         METADATA{:,i}.comps_num=EEG.Removed_ICA_Comps.comps_num;
    %           end
    %           if  isfield(EEG.Removed_ICA_Comps,'remcomp_ICA2')
    %         METADATA{:,i}.remcomp_ICA2=EEG.Removed_ICA_Comps.remcomp_ICA2;
    %           end
    %         end
    
    
%[EEG.event.urevent]=deal(EEG.event.epoch);
%EEG.urevent=EEG.urevent(strcmp({EEG.urevent.type},'128'));
%EEG.urevent=EEG.event; %EEG.urevent(1:size(EEG.event,2));
%EEG.urevent=EEG.urevent(1:size(EEG.event,2));
%EEG.epoch=EEG.METADATA{1, 7}.epochs{1, 7}
%METADATA{1,1}.epochs
% [EEG.epoch.event] = deal(EEG.event.urevent);
% [EEG.epoch.eventurevent] = deal(EEG.event.urevent);
% [EEG.epoch.eventtype] = deal(EEG.event.type);
% [EEG.epoch.eventlatency] = deal(EEG.event.latency);
% if size(EEG.epoch(1).eventlatency,2)>1
% for j=1:size(EEG.epoch,2)
% EEG.epoch(j).eventlatency=EEG.epoch(j).eventlatency{end};
% end
% end




%EEG = pop_editeventvals(EEG,'delete',[find(strcmpi({EEG.event.type},'old'))]);
%
% EEG.event = EEG.event(1,1);
% for i = 2:size(fileNames,1);
%     EEG.epoch = [EEG.event EEG.event(1,1)];
% end;
%
%
% EEG.epoch = EEG.epoch(1,1);
% for i = 2:size(fileNames,1);
%     EEG.epoch = [EEG.epoch EEG.epoch(1,1)];
% end;



%       EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp'); % Uri Alyagon 24032013 for lior

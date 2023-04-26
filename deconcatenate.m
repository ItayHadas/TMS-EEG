function deconcatenate
% need to have the EEG.Concatinated in the EEGset

%load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\chanlocs62.mat');
% clear EEG ALLEEG output fileNames pathName

[fileNames, pathName]=getANTFileNames;

for i=1:size(fileNames,1)
    if size(fileNames, 1)==1
        fileName=fileNames{i,1}';
    else
        fileName=fileNames{i,1};
    end;
    
[ num2str(i) ' out of ' num2str(size(fileNames,1)) '   -   ' fileName ]
% if isempty(EEG.Concatinated)
EEGcon = pop_loadset( [pathName fileName]);
%EEGcon = pop_loadset( [pathName fileNames{:}']);
EEGcon = eeg_checkset( EEGcon );

loc=[0 cumsum([EEGcon.Concatinated.EPOCHsize{:}])];

EEG=EEGcon; EEG.data=[]; EEG.setname=''; EEG = rmfield(EEG,'Concatinated')
for i=1:size(EEGcon.Concatinated.EPOCHsize,1)
    
    EEG.data=EEGcon.data(:,:,loc(i)+1:loc(i+1));
    EEG.setname=[ EEGcon.Concatinated.Filesnames{i}(:,1:end-4) '_DEcon'];
    EEG.filename = [EEG.setname '.set'];
    EEG.datfile = [EEG.setname '.fdt'];
    
    EEG.trials=EEGcon.Concatinated.trials{i};
    EEG.event=EEGcon.Concatinated.event{i};
    EEG.urevent=EEGcon.Concatinated.urevent{i};
    EEG.epoch=EEGcon.Concatinated.epoch{i};
    EEG.EPOCHsize=EEGcon.Concatinated.EPOCHsize{i};
    EEG.subject=EEGcon.Concatinated.subject{i};
    EEG.chanlocs=EEGcon.Concatinated.chanlocs{i};
    EEG.eventdescription=EEGcon.Concatinated.eventdescription{i};
    EEG.reject=EEGcon.Concatinated.reject{i};
    EEG.tmscut=EEGcon.Concatinated.tmscut{i};
    EEG.Removed_ICA_Comps=EEGcon.Concatinated.Removed_ICA_Comps{i};  
    EEG.Removed_ICA_Comps_while_Concatinated=EEGcon.Removed_ICA_Comps;	
if isfield(EEGcon.Removed_ICA_Comps,'Removed_Comps2')==1	

    EEG.Removed_ICA_Comps.concatenated_components=EEGcon.Removed_ICA_Comps.Removed_Comps2;
	end
	if isfield(EEGcon.Concatinated,'condition')==1
    EEG.condition=EEGcon.Concatinated.condition{i};
    end
    % concatenated ICA data removal
    EEG.icaact=[];  EEG.icawinv=[];  EEG.icasphere=[];  EEG.icaweights=[];  EEG.icachansind=[];
    
    
    %%
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, [pathName EEG.filename]);
end


end

end


clear all
[fileNames, pathName]=Z_getSetsFileNames('cdt');
display([pathName, fileNames{:}])
for i=1: size(fileNames,1)
    fileName=fileNames{i};
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
    EEG.session=EEG.setname(end-2:end);
    disp(['remove un-needed channels  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    remove_channels= {'lz' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
    EEG = pop_select( EEG,'nochannel',remove_channels);
    EEG = eeg_checkset( EEG );
    disp(['Downsampling 1000 Hz  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    EEG.data=double(EEG.data);
    EEG = pop_resample(EEG, 1000);
    disp(['Saving  - dataset ' num2str(i) '/' num2str(size(fileNames,1))]);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath', pathName ,'check', 'on','savemode','onefile','version','7.3');
end
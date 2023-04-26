function errrorsss = concatenate(pathName, subloc , condloc)
% needs the same number of electrodes per EEG.data
%load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\chanlocs62.mat');
% clear EEG ALLEEG output fileNames pathName


%   subloc=[2:3]; condloc=[5];

%pathName2=[pathName 'precont\']; %loading files from the last folder with the filtered files
pathName2=[pathName];
clear Dirfiles fileNames
Dirfiles=dir([pathName2 '*.set']);
%[fileNames, pathName]=getANTFileNames;
fileNames={Dirfiles.name}';
list1=char(fileNames); % list1(:,2:10)
tt=string(fileNames); %strfind
subs=unique(string(list1(:,subloc)));
errrorsss={};
for s=1:size(subs,1)
    subl=tt(find(strcmp(list1(:,subloc),subs(s))));
    subll=char(subl);
    %cond=unique(string(subll(:,condloc)));
    for p=str2num(unique(string(subll(:,condloc))))' % 1:3 %3 phases %1:size(unique(string(subll(:,condloc))),1) %
        cond1=subl(find(strcmp(string(subll(:,condloc)),num2str(p))));
        %cond1=cond(p);
        if size(cond1,1)==2
            
            
            ALLEEG = pop_loadset('filename',{cond1{:}},'filepath',strrep(pathName,'\','\\'));
            eponum=cellfun(@(x) size(x,3),{ALLEEG.data},'UniformOutput',false);
            EEG = pop_mergeset(ALLEEG,[1:size(eponum,2)]);
            EEG.Concatinated.EPOCHsize=eponum';
            EEG.Concatinated.Filesnames=cond1;
			EEG.Concatinated.condition={ALLEEG.condition};
            EEG.Concatinated.trials={ALLEEG.trials};
            EEG.Concatinated.event={ALLEEG.event};
            EEG.Concatinated.urevent={ALLEEG.urevent};
            EEG.Concatinated.epoch={ALLEEG.epoch};
            EEG.Concatinated.subject={ALLEEG.subject};
            EEG.Concatinated.chanlocs={ALLEEG.chanlocs};
            EEG.Concatinated.eventdescription={ALLEEG.eventdescription};
            EEG.Concatinated.reject={ALLEEG.reject};
            EEG.Concatinated.tmscut={ALLEEG.tmscut};
            EEG.Concatinated.Removed_ICA_Comps={ALLEEG.Removed_ICA_Comps};
            nam1=char(cond1(1));
            EEG.setname = ['S' nam1(2:5) '_concatenated_PRE_POST_singlepulse'];
            EEG.filename = [EEG.setname '.set'];
            EEG.datfile = [EEG.setname '.fdt'];
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, [pathName EEG.filename]);
            %out1=struct2table(rmfield(EEG.Concatinated,{'trials','event','urevent','epoch'}));
            clear -vars ALLEEG EEG eponum nam1
            
            
        elseif size(cond1,1)~=2
            'duplicate or missing file for the condition'
            cond1
            errrorsss(size(errrorsss,1)+1)={cond1};
            
        end
    end
end

errrorsss=cell2table(errrorsss');
writetable(errrorsss,strjoin([pathName 'errrorsss_concatinations.xlsx'])) ;






end




% cellstr(strfind(tt,subs(1)))
% isemty(strfind(tt,subs(1)))
% find(cell2mat(strfind(tt,subs(1))))

% ALLEEG = pop_loadset('filename',{fileNames{:}},'filepath',strrep(pathName,'\','\\'));
% eponum=cellfun(@(x) size(x,3),{ALLEEG.data},'UniformOutput',false);
% EEG = pop_mergeset(ALLEEG,[1:size(eponum,2)]);
% EEG.Concatinated.EPOCHsize=eponum';
% EEG.Concatinated.Filesnames=fileNames;
% EEG.Concatinated.trials={ALLEEG.trials};
% EEG.Concatinated.event={ALLEEG.event};
% EEG.Concatinated.urevent={ALLEEG.urevent};
% EEG.Concatinated.epoch={ALLEEG.epoch};
% EEG.setname = ['concatenated_conditions'];
% EEG.filename = [EEG.setname '.set'];
% EEG.datfile = [EEG.setname '.fdt'];
% EEG = eeg_checkset( EEG );
% EEG = pop_saveset( EEG, [pathName EEG.filename]);


% output=struct2table(rmfield(EEG.Concatinated,{'trials','event','urevent','epoch'}));
% writetable(output,[pathName 'output_concatenate_datasets.xlsx'])



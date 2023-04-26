function [ homoNum homoname] = homoElec( EEG, elecname )
%  converts electrodes namae to its number in spesific EEG dataset
%   elecname={ 'f3' 'f2' 'oz' 'fpz' 'f1' 'fc4' 'fc6' 'f4' 'fc4' 'o6' 'f5' 'fc3' 'fc5'  }
% elecname={  'fc4'  }
%  i=1;  EEG=ALLEEG(i)

%e = struct2cell(EEG.chanlocs);
e=struct2table(EEG.chanlocs);
exx=~(strcmpi('X',e.Properties.VariableNames) | strcmpi('Z',e.Properties.VariableNames));
e{:,[strcmp('double',varfun(@class,e,'OutputFormat','cell')) & exx]}=floor(abs(e{:,[strcmp('double',varfun(@class,e,'OutputFormat','cell')) & exx]}));
e{:,strcmp('double',varfun(@class,e,'OutputFormat','cell'))}=zscore(e{:,strcmp('double',varfun(@class,e,'OutputFormat','cell'))},0,1);
elect=table2cell(e(:,[true strcmp('double',varfun(@class,e(:,2:end),'OutputFormat','cell'))]));
%elect = [squeeze(e(1,1,:)) squeeze(e(3,1,:)) squeeze(e(4,1,:)) squeeze(e(5,1,:)) squeeze(e(6,1,:)) squeeze(e(7,1,:)) squeeze(e(8,1,:)) squeeze(e(9,1,:))];
elect = elect(find(~cellfun('isempty', elect(:,2))),:);
locmat=cell2mat(elect(:,2:end));
locmat(:,std(locmat,0,1)==0)=[];
%locmat=cell2mat(e(:,2:end));
[neighbors, ~]=rangesearch(locmat,locmat,0.6);
%rangesearch(cell2mat(e(:,2:end)),cell2mat(e(:,2:end)),6)
% varfun(@class,e),'OutputFormat','logic')
eNum = [];
homoname = {};
homoNum = [];

for a = 1:size(elecname,2)
    if isempty(find(strcmpi(elecname{1,a},elect(:,1))))
        
        warning(['electrode ' elecname{1,a} ...
            ' is not present in the EEG.chanlocs set. Electrode not included ']);
        homoname(1,size(homoname,2)+1) = {'NONEabab'};
         homoNum(1,size(homoNum,2)+1) =  [nan];
         eNum (1,size(eNum,2)+1) =[nan]  ;
    elseif ~isempty(find(strcmpi(elecname{1,a},elecname(1:a-1))))
        ff=find(strcmpi(elecname{1,a},elect(:,1)));
        warning(['electrode ' EEG.chanlocs(ff).labels ' on channel ' num2str(EEG.chanlocs(ff).urchan)...
            ' is a duplicate. Electrode not included ']);
        homoname(1,size(homoname,2)+1) = {'NONEabab'};
         homoNum(1,size(homoNum,2)+1) =  [nan];
         eNum (1,size(eNum,2)+1) =[nan]  ;
    else
        hom=[];
        eNum (1,size(eNum,2)+1) =  find(strcmpi(elecname{1,a},elect(:,1)));
        hom=cell2mat(neighbors(eNum(1,a)));%  hom=cell2mat(neighbors(43));
        %hom=hom(hom~=eNum);
        % elecname
        if size(hom,2)>1
            
            homoname(1,size(homoname,2)+1)= {EEG.chanlocs(hom(hom~=eNum(end))).labels};   %hom(strcmpi(hom,elecname(a))==0);
            if ~isempty(EEG.chanlocs(hom(hom~=eNum(end))).urchan)
            homoNum(1,size(homoNum,2)+1)= EEG.chanlocs(hom(hom~=eNum(end))).urchan  ; %find(strcmpi(homoname{1,a},elecc(:,1)));
            else
             homoNum(1,size(homoNum,2)+1)= hom(hom~=eNum(end)) ;   
            end
        elseif    size(hom,2)==1
            
            homoname(1,size(homoname,2)+1) = {'NONEabab'}; %{EEG.chanlocs(hom(hom~=eNum)).labels} ;%elecc(eNum(1,a),1);
            homoNum(1,size(homoNum,2)+1) =  [nan]; %EEG.chanlocs(hom(hom~=eNum)).urchan ; %eNum(1,a);
            warning (['electrode ' EEG.chanlocs(eNum(a)).labels ' on channel ' num2str(EEG.chanlocs(eNum(a)).urchan)  ' has no Homologue']) % elecc(eNum(1,a),1)  - but will remain in the vector
            
        end
    end
    % clear eNum eLoc hom homoname homoNum
end
homoname=homoname(~strcmp(homoname,'NONEabab'));
homoNum=homoNum(~isnan(homoNum));





end


%% old


%
%elecc(:,1)=elect(:,1) ;
%elecc(:,2)=cellstr(locmat);



%elect = [squeeze(e(1,1,:)) squeeze(e(3,1,:)) squeeze(e(5,1,:)) squeeze(e(6,1,:)) squeeze(e(7,1,:))];

% locmat = num2str(floor(abs(cell2mat(elect(:,2:7)))));
% floor(abs(cell2mat(elect(:,2:7))))
% num2str(floor(cell2mat(elect(:,2:7))))
% sum(floor(cell2mat(elect(:,2:7))),2)
% jj=sum(floor(abs(cell2mat(elect(:,2:7)))),2)
% jj=floor(abs(cell2mat(elect(:,2:7))))
% str2mat(strrep(string(num2str(floor(abs(cell2mat(elect(:,2:5)))))),' %','')); %old version
% clear -vars elec


%         eLoc=elecc(eNum(1,a),2);
%         hom=elecc(strcmp(elecc(:,1),eLoc),1);


% display (strjoin(['electrode ' EEG.chanlocs(hom(1)).labels ' - ' num2str(EEG.chanlocs(hom(1)).urchan)  ' has no Homologue - but will remain in the vector'])) % elecc(eNum(1,a),1)




function v = LOCF(locq)
%locq=locd
%locq=locn
for j=1:size(locq,1)
    %class(locq)
    %locq=num2cell(locq)
if  strcmp(class(locq{j,:}),'double') %strcmp(class(locq),'double') %
    idx = (~isnan(locq{j,:})); %non nans
    vr{j,:}=locq{j,:}(:,idx);
    if ~isempty(vr{j,:})
        vr{j,:}=[nan vr{j,:}];
        v{j,:}=vr{j,:}(cumsum(idx)+1);
        else
        v{j,:}=locq{j,:};
    end
    
elseif strcmp(class(locq{j,:}),'cell') %strcmp(class(locq),'cell') %
    
    locq{j,:}(cell2mat( cellfun(@(x)any(isnan(x)),locq{j,:},'UniformOutput', false)))={''};
    idx =cellfun(@isempty, locq{j,:}  )==0 ;
    vr{j,:}=locq{j,:}(:,idx);
    if ~isempty(vr{j,:})
        vr{j,:}=[{''} vr{j,:}] ;
        v{j,:}=vr{j,:}(cumsum(idx)+1);
        else
        v{j,:}=locq{j,:};
    end
   
    elseif strcmp(class(locq{j,:}),'string') %strcmp(class(locq),'cell') % 
    
    %locq{j,:}(cell2mat( cellfun(@(x)any(isnan(x)),locq{j,:},'UniformOutput', false)))={''};
    idx =cellfun(@isempty, locq{j,:}  )==0 ;
    vr{j,:}=locq{j,:}(:,idx);
    if ~isempty(vr{j,:})
        vr{j,:}=[{''} vr{j,:}] ;
        v{j,:}=vr{j,:}(cumsum(idx)+1);
        else
        v{j,:}=locq{j,:};
    end
end



%for j=1:size(locq,1)
    % clear vr v
%    vr{j,:}=locq(j,idx(j,:));
%     if ~isempty(vr{j,:})
%         if strcmp(class(locq{j,:}),'double')
%             vr{j,:}=[nan vr{j,:}];
%         elseif strcmp(class(locq{j,:}),'cell')
%             vr{j,:}=[{''} vr{j,:}] ;
%         end
%         v{j,:}=vr{j,:}(cumsum(idx)+1);
%     else
%         v{j,:}=locq{j,:}
%     end
% end

end
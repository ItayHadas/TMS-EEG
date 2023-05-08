%dropOFFchan
function [DROP,  dropNaN, dropdead, drophigh] = dropOFFchan( EEG )
%drophigh dropdead

timeind = EEG.times > -300 & EEG.times < -20 | EEG.times> 20 & EEG.times < 350; % time of interest
%drophigh=nan;
dropNaN=nan;

 % pop_eegplot( EEG, 1, 1, 1);
for n = 1:size(EEG.data,1)
    if sum(sum(isnan(EEG.data(n,:,:))))>0
       dropNaN=[dropNaN n];
    end
end
%%
for c = 1:size(EEG.data,1)
    vec=reshape(EEG.data(c,timeind,:),[size(EEG.data(c,timeind,:),2)*size(EEG.data(c,:,:),3) 1]);
    vec2=reshape(EEG.data(c,:,:),[size(EEG.data(c,:,:),2)*size(EEG.data(c,:,:),3) 1]);
    %el(c)=std(abs(vec));
    %elm(c)=mean(abs(vec));
    elmed(c)=median(abs(vec));
    elmed2(c)=median(abs(vec2));
    
    %clear -vars vec
end


dropdead=find(elmed2<0.01);

eyes=elecName(EEG,{'fp1' 'fpz' 'fp2'});
noeyes=ones(1,size(EEG.chanlocs,2)); noeyes(1,eyes)=0; noeyes=logical(noeyes); a=1:size(EEG.data,1);a(noeyes);
drophigh=find(zscore(elmed)>5);
drophigh=unique([drophigh find(zscore(elmed2)>9)]);
drophigh=intersect(a(noeyes),drophigh);
%drophigh= [drophigh find(zscore(elmed)>5)];

% figure; plot(el); title('STD');
% figure; plot(elm); title('Mean');
% figure; plot(elmed); title('Median');
% 
% figure; hist(zscore(el)); title('STD');
% figure; hist(zscore(elm)); title('Mean');
% figure; hist(zscore(elmed)); title('Median');


%%


%dropNaN=dropNaN(2:end);
%dropdead=[] %canceling the drop silent/DEAD channels
% drophigh=drophigh(2:end);
%DROP=[ dropNaN drophigh dropdead ]; %drophigh 
DROP=[  drophigh dropdead ]
%EEG=pop_select(EEG,'nochannel', DROP);
%clear c d DROP elmed


%% OLD
% %if mean(std((EEG.data(c,:,:))))<0.1
    %   dropdead=[dropdead c];
    %end
% for d = 1:size(EEG.data,1)
%     if mean(std((EEG.data(d,:,:))))>mean(mean(std(EEG.data)))*7
%        drophigh=[drophigh d];
%     end
% end

function comprem(EEG,channels,components,comp)
% look at the TEP before and after specifyed components removal
%   channels=[7 8];%cahnnels to look at
%   components=[1 2 3 4 6 11 5 7 8 9 10];%components to remove [1 2 3 4 6 11]

if isempty(channels)
    error 'channels vector is empty';
    return
end
if isempty(components)
    error 'components vector is empty';
    return
end

EEG1=pop_subcomp(EEG,components,0);
%if strcmp(chanplot,'on')
a=figure;  subplot(2,1,1); hold on;
title(['avearege  ',num2str(channels),' after removing components:',num2str(components)])

plot(EEG.times,[permute(mean(EEG.data(channels,:,:),1),[2 3 1])],'color',[ 0 0 0 0.1])
plot(EEG1.times,permute(mean(EEG1.data(channels,:,:),1),[2 3 1]),'color',[ 1 0 0 0.4])
xlim([-200 400]); ylim([-50 50]); ax=gca; ax.XTick = [ -100 0 30 60 100 200 250 300];
subplot(2,1,2);hold on;
plot(EEG.times, squeeze(mean(mean(EEG.data(channels,:,:),3),1)),'color',[ 0 0 0 1])
plot(EEG1.times, squeeze(mean(mean(EEG1.data(channels,:,:),3),1)),'color',[ 1 0 0 1])
xlim([-200 400]); ylim([-20 20]);ax=gca; ax.YAxisLocation='origin'; ax.XAxisLocation='origin'; ax.XTick = [ -100 0 30 60 100 200 250 300];
legend({['Before'], ['After']})
hold off;

%
%
if strcmp(comp,'on')
    ALLEEG1(1)=EEG; ALLEEG1(1).setname='original';
    ALLEEG1(2)=EEG1; ALLEEG1(2).setname=['components removed: ' num2str(components)];
    %pop_comperp( ALLEEG1, 1, [1 2] ,[],'tlim',[-200 400], 'addavg','off','addstd','off','addall','on','diffavg','off','diffstd','off','tplotopt',{'ydir' 1});
    pop_comperp( ALLEEG1, 1, [1 2],[],'tlim',[25 350],'chans',[ 1:size(EEG.chanlocs,2)],'addavg','off','addstd','off','addall','on','diffavg','off','diffstd','off','tplotopt',{'ydir' 1});
%[ 53 39 1 19 33 5 47 13 55 27 41 23 2 15 49 7 35 20]
end
end


function [EEG] = FFT_comp(EEG, minfreq, maxfreq )
% minfreq=59; maxfreq=62;


freq=2:0.5:130;
freqind=freq>=minfreq & freq<=maxfreq;
%freqind2=(freq>=minfreq-20 & freq<minfreq) | (freq>maxfreq & freq<=maxfreq+20);
freqind2=(freq>=1 & freq<50); % | (freq>maxfreq & freq<=maxfreq+20);

  
%% calculating FFT and removing the component with prominent minfreq-maxfreq window

%Ranks components, outputs variance
[EEG, ~] = tesa_sortcomps(EEG); 
% pop_plotdata(EEG, 0, [1:size(EEG.icaweights,1)] , [1:size(EEG.data,3)], EEG.setname, 0, 1, [0 0]);
%Calculates time course and variance
%compTimeCourse = arrayfun(@(x)eeg_getdatact(EEG, 'component', [x], 'projchan', []),1:size(EEG.icawinv,2),'UniformOutput', false);
comps = size(EEG.icaweights,1);
%comps=10
if comps>30; comps=30; end
compTimeCourse = arrayfun(@(x)eeg_getdatact(EEG, 'component', [x], 'projchan', [1:size(EEG.data,1)]),[1:comps],'UniformOutput', false);

%%
%comps=10
if comps>30; comps=30; end
for compNum =1:comps
    %Calculates time course
    [~, elec]=max(abs(EEG.icawinv(:,compNum))); 
    temp = compTimeCourse{1,compNum};
    temp = squeeze(temp(elec,:,:));
    %Calculates FFT (uV/Hz)
    T = 1/EEG.srate;             % Sample time 
    L = size(EEG.times,2);       % Length of signal
    NFFT = 2^nextpow2(L);        % Next power of 2 from length of y
    f = EEG.srate/2*linspace(0,1,NFFT/2+1); % Frequencies
    y = reshape(temp,1,[]);
    Y = fft(zscore(y),NFFT)/L;
    Yout = (abs(Y).^2);
    %[~, index]=min(abs(f-(0.5)));
    
    for x = 1:size(Yout,1)
        for a=1:size(freq,2)
            [~, index1]=min(abs(f-((freq(1,a)-0.25))));
            [~, index2]=min(abs(f-((freq(1,a)+0.25))));
            Y2(x,a)=mean(Yout(x,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
        clear -vars index1 index 2 ttt
        end
    end
    %
    
    %[val1,winF1] = min(abs(minfreq)); % removed on 13/05/2019 by itay
    %[val2,winF2] = min(abs(maxfreq)); % removed on 13/05/2019 by itay
    
    %winFreq = mean(Y2(:,winF1:winF2),2);
    %winFreq = mean(Y2(:,floor(val1/2):ceil(val2/2)),2);
    %winFreq = max(Y2(:,floor(val1/2):ceil(val2/2)));
    winFreq = max(Y2(:,freqind));
    %allFreq = mean(Y2(:,[[floor(val1/2)-10]:[floor(val1/2)] [ceil(val2/2)]:[ceil(val2/2)+10]]),2);
    allFreq = mean(Y2(:,freqind2),2);
    Ratio(compNum) = winFreq./allFreq;
    %Ratio(compNum) = winFreq-allFreq;
%     figure;  plot(freq,Y2);
%     title(['component no. ' num2str(compNum) '  ratio:' num2str(Ratio(compNum))]);;

%clear -vars val1 winF1 val2 winF2 winF2 allFreq temp T L NFFT f y Y Yout c index index1 index2 Y2 a x
clear -vars allFreq temp T L NFFT f y Y Yout c index Y2 a x
end
% figure; hist(zscore(Ratio))

    %%
%remcomps=find(zscore(Ratio)>=2.5)%Ratio(find(zscore(Ratio)>=2));
disp('Line noise components removed:')
remcomps=find(Ratio>=200)
%
% Ratio(find(Ratio>=100))
% TMS_pop_selectcomps(EEG, [remcomps] );
% pop_eegplot( EEG, 0, 1, 1, {}, 'spacing', 40, 'dispchans', 10);
%%
% pop_plotdata(EEG, 0, [1:size(EEG.icaweights,1)] , [1:size(EEG.data,3)], EEG.setname, 0, 1, [0 0]);
if sum(remcomps)>0
    
    setname1=EEG.setname;
    filename1=EEG.filename;
    datfile1=EEG.datfile;
    EEG=pop_subcomp( EEG, [remcomps], 0); % remove the freq window components
    EEG.setname=setname1;
    EEG.filename=filename1;
    EEG.datfile=datfile1;
    %EEG.setname=setname;
    EEG = Z_append(EEG,'_freq');
    EEG.Removed_ICA_Comps.FFTcomp={num2str(remcomps)};
else
    EEG = Z_append(EEG,'_NOfreq');
end

EEG = eeg_checkset( EEG );
%EEG = pop_saveset( EEG, [pathName 'ICA1\' EEG.filename]);

end

    

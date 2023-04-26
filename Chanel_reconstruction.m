function draw_specto (comp, electrodes, timewin, freqwin,powerlim)


%%
% electrodes= {'fc3' 'fc5'}
% powerlim=[0 2.5]
% timewin=[-400 1000];
% freqwin=[0 80];
% comp='ersp' %'itc'
[fileNames,pathName] = uigetfile({'*.mat'},'Choose Wavelet File','on');

load([pathName, fileNames]);


timeind=grandWave(1,1).times>=timewin(1) & grandWave(1,1).times<=timewin(2);
freqind=grandWave(1,1).freqs>=freqwin(1) & grandWave(1,1).freqs<=freqwin(2);

for j=1:size(electrodes,2)
    elec(j)=find(strcmpi(electrodes(j),{grandWave.wavelet.channelLabel}));
end

for jj=1:size(elec,2)
    if strcmp(comp,'itc')
        ersp(:,:,jj)=mean([grandWave.wavelet([elec(jj)]).itc], 3);
    elseif strcmp(comp,'ersp')
        ersp(:,:,jj)=mean([grandWave.wavelet([elec(jj)]).ersp], 3);
        %     elseif strcmp(comp,'blcpow')
        %         ersp(:,:,jj)=mean([grandWave.wavelet([elec(jj)]).powbase], 3);
    end
    
end
ersp=mean(ersp,3);


figure;
imagesc(grandWave(1,1).times(timeind), grandWave(1,1).freqs(freqind), ersp(freqind,timeind),powerlim);
if strcmp(comp,'itc')
    title(strrep(['ITC' fileNames],'_',' '))
elseif strcmp(comp,'ersp')
    title(strrep(['ERSP' fileNames],'_',' '))
end
colorbar; colormap jet;
set(gca, 'YDir','normal');

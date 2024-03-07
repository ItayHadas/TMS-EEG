% matlab -nosplash -nodesktop -r "workspace"
 
clear all
if (ispc)
    sep='\';
    not_sep='/';
    rep_space = ' ';
    GITS='D:\GITs\';
    Path = dir('A:\WorkingSet\WellcomeLeap_TEP\**\*SPD_*.set');
    %Path = dir('A:\WorkingSet\TEPs_GRANDS')
    outdir='A:\WorkingSet\TEPs_GRANDS';
    eeglabdir='D:\MATLAB\EEGLAB'; 
    FTdir='D:\MATLAB\fieldtrip'
elseif (ismac || isunix)
    sep='/';
    not_sep='\';
    rep_space = '\ ';
    GITS='/media/ipp/DATA/GITs';
    Path = dir('/mnt/INTERPSYC/DATA/Bipolar_800601/Neurophysiology_Data/**/*SPD_*.cdt');
    outdir='/media/ipp/DATA/EEG_DATA/Bipolar_TEP'; 
    mkdir outdir
    eeglabdir='/media/ipp/DATA/Documents/MATLAB/eeglab2023.1'; 
    FTdir='/media/ipp/DATA/Documents/MATLAB/'
    addpath('/media/ipp/DATA/GITs/Localization/','/media/ipp/DATA/GITs/TMS-EEG/' )
end
run([FTdir sep 'ft_defaults'])
addpath([GITS sep 'TESA']); 
addpath(genpath([GITS sep 'TMS-EEG']),genpath([GITS sep 'AARATEPPipeline']));
addpath(eeglabdir); 
eeglab



EEG{1}=pop_loadset('filepath' , Path(7).folder, 'filename',Path(7).name);

color=brewermap(2,'set1');
%figure; topoplot([],EEG{2}.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG{2}.chaninfo)
    
%%
elec={'fc3' 'fc4' 'f3' 'f5' 'f1' 'fc1' 'c1' 'c3' 'c5'}
%elec={ 'fc3' 'fc4' 'f3' 'f5'}
timel=[-900 900];
LMFP=[]; LMFPste=[]; LMFPM=[];
for ii=1:length(EEG)
timeInd=EEG{ii}.times>=timel(1) & EEG{ii}.times<=timel(2);
[elecnum, elecname]=elecName(EEG{ii},elec);
LMFP(:,:,ii)=double(squeeze(std(EEG{ii}.data(elecnum,timeInd,:),0,1)));
end
[Y,E] = discretize(LMFP(1,:),size(LMFP,2))
LMFPM=smoothdata(mean(LMFP(:,[1:30]),2,"omitnan"),'gaussian',40);
figure; plot(EEG{1}.times(timeInd),LMFPM(:,:),'LineWidth',2); hold on;
hline(mean(LMFPM(:,:)));

figure; tiledlayout(1,10)
for tt=1:10
LMFPM=smoothdata(mean(LMFP(:,[1:tt+2]),2,"omitnan"),'gaussian',40);
 nexttile; plot(EEG{1}.times(timeInd),LMFPM(:,:),'LineWidth',2);
end

%%

figfig1=figure('position', [500 200 900 400]); hold on;
pp1=patch([EEG{1}.times(timeInd) fliplr(EEG{1}.times(timeInd))],[LMFP(:,:,1)-LMFPste(:,:,1) fliplr(LMFP(:,:,1)+LMFPste(:,:,1))],color(1,:),'EdgeColor',[color(1,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %
p1=plot(EEG{1}.times(timeInd),LMFP(:,:,1),'LineWidth',2);
p1.Color=color(1,:);
pp2=patch([EEG{2}.times(timeInd) fliplr(EEG{2}.times(timeInd))],[LMFP(:,:,2)-LMFPste(:,:,2) fliplr(LMFP(:,:,2)+LMFPste(:,:,2))],color(2,:),'EdgeColor',[color(2,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
p2=plot(EEG{2}.times(timeInd),LMFP(:,:,2),'LineWidth',2);
p2.Color=color(2,:);
yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
%f=fill([statwin statwin(2) statwin(1)],[yy yy 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
chH = get(gca,'Children');
set(gca,'Children',[chH(2:end) ;chH(1)])
%fb=fill([mintrim maxtrim maxtrim mintrim],[yy yy 0 0],'w','edgecolor','none');
%axis([timeWin 0 amp]);
ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
legend('show'); clear title; 
% title(['Local Mean Field Power ' strjoin([elecname]) '  time:'...
%     strrep(num2str(statwin),'  ','_')], 'interpreter', 'none'); 
title(['Local Mean Field Power ' strjoin([elecname]) ...
    ], 'interpreter', 'none'); 
ylabel ('\muv'); xlabel ('milliseconds');
legend([p1 p2],{strrep([EEG{1}.setname ' n=' num2str(size(EEG{1}.subjects,2))],'_',' '),...
strrep([EEG{2}.setname ' n=' num2str(size(EEG{2}.subjects,2))],'_',' ')})
ylim([min(get(gca,'YLim')) round(max([pp1.YData ; pp2.YData]))+1])
set(ax ,'Layer', 'Top')
% ['ADHD  n=' num2str(size(corrsADHD(NnA,:),1))]
hold off;
clearvars f p1 p2 ax

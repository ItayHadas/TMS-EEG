function [ stats head TEPVec1 g] = TEPstats(ALLEEG, dataset, Groups, name, fixfit, facelift,rem_outlie , violinplt, ISP, rectyfied, wavgrph,LMFPgrph,GMFPgrph, bargrph, timeWin, statwin, stichsize, elec, amp, method, peaktype, smallWinHalfLength)
%% options optsTEP
%optsTEP
g={};
% headminmax
%clearvars -except fileNames pathName ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG LASTCOM PLUGINLIST STUDY eeglabUpdater Final

%% TEP plot
%ADHDcolor=[238/255 17/255 17/255];%'r'; %[0 0 0];
%Healthycolor=[43/255 87/255 151/255]; % 'b' ;%[166/255 166/255 166/255];
%condName='TEP';
% co=brewermap(3,'RdBu');
% ADHDcolor=co(1,:);%'r'; %[0 0 0];
% Healthycolor=co(3,:); % 'b' ;%[166/255 166/255 166/255];
co=brewermap(4,'PuOr'); %PuOr RdBu
ADHDcolor=co(1,:);
%ADHDcolor1=co(2,:);%'r'; %[0 0 0];
Healthycolor=co(4,:); % 'b' ;%[166/255 166/255 166/255];
%Healthycolor1=co(3,:);

Groups=[Groups {ADHDcolor Healthycolor}];

ISP_lag=10;

if statwin(1)-statwin(2)>=0
    warning('stat time window is BAD')
    return
end
if timeWin(1)-timeWin(2)>=0
    warning('Graph time window is BAD')
    return
end

%dat=[1 2]  dataset=[1 2]
color{dataset(1)}=ADHDcolor;
color{dataset(2)}=Healthycolor;

for i=dataset %1:size(ALLEEG,2)
    
    %  i=1
    if strfind(ALLEEG(i).filename,Groups{1})
        color{i}=Groups{3};
    elseif strfind(ALLEEG(i).filename,Groups{2})
        color{i}=Groups{4};
    end
    timeInd=ALLEEG(i).times>=timeWin(1) & ALLEEG(i).times<=timeWin(2);
    timeInd_fix=ALLEEG(i).times>=-400 & ALLEEG(i).times<=0;
    statInd=ALLEEG(i).times>=statwin(1) & ALLEEG(i).times<=statwin(2);
    [~, fixpre] = min(abs(ALLEEG(i).times-(-500)));
    fixIndpre=ALLEEG(i).times>=-500 & ALLEEG(i).times<=-1;
    [~, fixpost] = min(abs(ALLEEG(i).times-(0)));
    fixIndpost=ALLEEG(i).times>=16 & ALLEEG(i).times<=600;
    stich=ALLEEG(i).times>=-2 & ALLEEG(i).times<=stichsize+2;
    % elecnum=[]; elecname={};
    [elecnum, elecname]=elecName(ALLEEG(i),elec);
    %    [elecHOMOnum, elecHOMOname]=homoElec(ALLEEG(i),elec);
    % [elecHOMOnum, elecHOMOname]=homoElec(ALLEEG(i),'fc3')
    % smoothing the stich
    t=ALLEEG(i).data;
    %elecs=unique([elecnum elecHOMOnum]);
    %t=FilterTEP(t, [0.01], [30],elecs,ALLEEG(i).srate);
    %t(:,stich,:)=sgolayfilt(double(t(:,stich,:)),3,11);
    %   for sub=1:size(t,3)
    %    t(:,stich,sub)=hampel(double(t(:,stich,sub)));
    %end
    %t(:,stich,:)=medfilt1(double(t(:,stich,:)),2);
    %t(:,stich,:)=smoothdata(double(t(:,stich,:)));
    %
    if strcmp(facelift,'on')
        %smoo_method='movmedian'; smoo_winsize=4; % movmean movmedian
        %t(elecnum,timeInd_fix,:)=filloutliers(t(elecnum,timeInd_fix,:),'clip','movmean',5,'ThresholdFactor',1);
        %t(elecnum,timeInd_fix,:)=t(elecnum,timeInd_fix,:)./2;
        % TEPVec1(:,:,i)=detrend(TEPVec1(:,:,i),1);
        %TEPVec1(:,:,i)=double(nanmean(mean(t(elecnum,timeInd,:),1),3));
        
        % t=FilterTEP(t, [0.01], [30],elecnum,ALLEEG(i).srate);
        % [zl, pl, kl]=butter(50,40./(ALLEEG(i).srate/2)); %Lowpass
        % [x_filtlow,y_filtlow] = zp2sos(zl,pl,kl);
        %  for elecs=elecnum %1:size(EEG.data,1)
        %         for epo=1:size(t,3)
        %             %EEG.data(elec,:)=filtfilt(x_filt,y_filt,double(EEG.data(elec,:)));
        %             %K2(elec,:)=filtfilt(x_filt,y_filt,double(K2(elec,:)));
        %             %K2(elec,:,epo)=filtfilt(x_filthigh,y_filthigh,K2(elec,:,epo));
        %             t(elecs,:,epo)=filtfilt(x_filtlow,y_filtlow,t(elecs,:,epo));
        %            % K2(elec,:,epo)=filtfilt(x_filtstop,y_filtstop,K2(elec,:,epo));
        %         end
        %     end
        
        TEPVEC(:,:,i)=double(nanmean(mean(t(elecnum,:,:),1),3));
        
        FixVec1pre(:,:,i)=double(nanmean(mean(TEPVEC(:,fixIndpre,:),1),3))-TEPVEC(:,fixpre,i);
        FixVec1post(:,:,i)=double(nanmean(mean(TEPVEC(:,fixIndpost,:),1),3))-TEPVEC(:,fixpost,i);
        % Fit model to data. [fitresultpre, goodness, output, convmsg]
        
        %[0:size(FixVec1pre(:,:,i),2)-1]', FixVec1pre(:,:,i)'
        [xData, yData] = prepareCurveData( [0:size(FixVec1pre(:,:,i),2)-1]', FixVec1pre(:,:,i)' );
        
        % Set up fittype and options.
        ft = fittype( 'rat55' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Algorithm = 'Levenberg-Marquardt';
        opts.Display = 'Off';
        opts.MaxFunEvals = 500;
        opts.MaxIter = 500;
        opts.Normalize = 'on';
        opts.Robust = 'LAR';
        opts.StartPoint = [0 0.197278942097314 0.111184983788934 0.297354293368622 0.396418535561727 0.420755683348492 0.311475357650453 0.693843170932226 0.0918718332714995 0.402088619037046 0.295180804029962];
        opts.TolFun = 1e-10;
        opts.TolX = 1e-10;
        
        % Fit model to data.
        [fitresultpre, ~] = fit( xData, yData, ft, opts );
        
        
        %          fitresultpre=fit([0:size(FixVec1pre(:,:,i),2)-1]', FixVec1pre(:,:,i)', 'poly3');% ALLEEG(i).times(timeInd)'
        
        %polyfix([0:size(FixVec1pre(:,:,i),2)-1], FixVec1pre(:,:,i)',[],[],[])
        
        % p = polyfix(x,y,n,xfix,yfix)
        %fitresultpost=fit([0:size(FixVec1post(:,:,i),2)-1]', FixVec1post(:,:,i)', 'poly2' );
        %   figure; plot(feval(fitresult,ALLEEG(i).times(timeInd))')
        % figure;
        TEPVEC(:,:,i)=TEPVEC(:,:,i)-TEPVEC(:,fixpre,i);
        TEPVEC(:,fixIndpre,i)=[TEPVEC(:,fixIndpre,i)-feval(fitresultpre,[0:size(FixVec1pre(:,:,i),2)-1])'];
        %TEPVEC(:,:,i)=TEPVEC(:,:,i)-TEPVEC(:,fixpost,i);
        %%
        %TEPVEC(:,fixIndpost,i)=[TEPVEC(:,fixIndpost,i)-feval(fitresultpost,[0:size(FixVec1post(:,:,i),2)-1])'];
        samp1=15.*[round(ALLEEG(i).srate/1000)]; samp2=15 .*[round(ALLEEG(i).srate/1000)];
        xx=TEPVEC(:,[fixpost-samp1:fixpost+samp2],i);
        xx(5:end-5)=nan;
        inn = ~isnan(xx);
        i1 = (1:numel(xx)).';
        pp = interp1(i1(inn),xx(inn),'pchip','pp'); % spline 'linear'
        out_1 = fnval(pp,linspace(i1(1),i1(end),size(xx,2)));
        TEPVEC(:,[fixpost-samp1:fixpost+samp2],i)=out_1;
        
        %TEPVEC(:,:,i)=smoothdata(TEPVEC(:,:,i),'gaussian',10);
        TEPVec1(:,:,i)=double(smoothdata(TEPVEC(:,timeInd,i),'gaussian',7));
        
        
        TEPVecste(:,:,i)=double(ste(smoothdata(mean(t(elecnum,timeInd,:),1),'gaussian',7),3));
        % plot(TEPVEC(:,fixIndpre,i)), hold on
        % plot(feval(fitresultpre,[0:size(FixVec1pre(:,:,i),2)-1]'))
        % plot([TEPVEC(:,fixIndpre,i)-feval(fitresultpre,[0:size(FixVec1pre(:,:,i),2)-1])'])
        %    else
        %        TEPVEC(:,:,i)=double(nanmean(mean(t(elecnum,:,:),1),3));
        %        %[~,rate] = rat(1000/ALLEEG(1).srate, 1e-12);
        % [~, ind_t0] = min(abs(ALLEEG(i).times-(0)));
        
        
        %TEPVEC(:,:,i)=double(nanmean(smoothdata(mean(t(elecnum,:,:),1),'gaussian',7),3)); %movmean gaussian ,rate
        
        
        
        % TEPVec1(:,:,i)=[TEPVec1(:,:,i)-feval(fitresult,ALLEEG(i).times(timeInd))'];
        %         TEPVec1(:,:,i)=double(nanmean(mean(t(elecnum,timeInd,:),1),3));
    else
        
        % TEPVec1(:,:,i)=double(TEPVEC(:,timeInd,i));
        TEPVec1(:,:,i)=double(mean(mean(t(elecnum,timeInd,:),1),3));
        TEPVecste(:,:,i)=double(ste(mean(t(elecnum,timeInd,:),1),3));
    end
    
    % TEPVecste(:,:,i)=double(ste(smoothdata(mean(t(elecnum,timeInd,:),1),'gaussian',7),3));
    %smoothdata(mean(t(elecnum,timeInd,:),1),50,'movmean',20)
    
    %
    if strcmp(fixfit,'on')
        TEPVec1(:,:,i)=detrend(TEPVec1(:,:,i),1);
        % Fit model to data.
        fitresult=fit([0:size(TEPVec1(:,:,i),2)-1]', TEPVec1(:,:,i)', 'poly2'  );% ALLEEG(i).times(timeInd)'
        %   figure; plot(feval(fitresult,ALLEEG(i).times(timeInd))')
        % figure;
        TEPVec1(:,:,i)=[TEPVec1(:,:,i)-feval(fitresult,ALLEEG(i).times(timeInd))'];
        % figure; plot(TEPVec1(:,:,i))
    end
    
    
    
    
    %% computations - GMFP, LMFP, ISP, rectyfied
    %TEPstd(:,:,i)=
    
    if strcmp(GMFPgrph,'on')
        % gmpf(i,j,k,ALLEEG(ind).subjects(p)-100,q) =sqrt(sum((temp(:,q,p)-mean(temp(:,q,p),1)).^2)/63);%Global Mean Field Power - sam's quation
        %GMFP(:,:,i)=double(mean(sqrt(sum((t(:,timeInd,:)-mean(t(:,timeInd,:),1)).^2)/size(t,1)),3)); %Global Mean Field Power
        %elecnum_G=size(ALLEEG(i).data,1);
        
        GMFP(:,:,i)=double(squeeze(nanmean(std(t(:,timeInd,:),0,1),3)));
        %double(squeeze(nanmean(sqrt(nansum([t([1:elecnum_G],timeInd,:)-nanmean(t([1:elecnum_G],timeInd,:),2)].^2,1)./size([1:elecnum_G],2)),3)));
        GMFP(:,:,i)=smoothdata(GMFP(:,:,i),'gaussian',7);
        GMFPste(:,:,i)=double(squeeze(ste(std(t(:,timeInd,:),0,1),3)));
        %double(squeeze(ste(sqrt(nansum([t([1:elecnum_G],timeInd,:)-nanmean(t([1:elecnum_G],timeInd,:),2)].^2,1)./size([1:elecnum_G],2)),3)));
        GMFPste(:,:,i)=smoothdata(GMFPste(:,:,i),'gaussian',7);
        %GMFPste(:,:,i)=double(ste(sqrt(sum((t(:,timeInd,:)-mean(t(:,timeInd,:),1)).^2)/size(t,1)),3));
    end
    
    if strcmp(LMFPgrph,'on') %& size(elecnum,2)>1
        % gmpf(i,j,k,ALLEEG(ind).subjects(p)-100,q) =sqrt(sum((temp(:,q,p)-mean(temp(:,q,p),1)).^2)/63);%Global Mean Field Power - sam's quation
        % subsT.stat=  double(squeeze(nanmean(sqrt(nansum([ALLEEG(i).data(elecnum,statInd,:)-nanmean(ALLEEG(i).data(elecnum,statInd,:),2)].^2,1)./size(elecnum,2)),2)));
        LMFP(:,:,i)=double(squeeze(nanmean(std(t(elecnum,timeInd,:),0,1),3)));
        % double(squeeze(std(nanmean(t(elecnum,timeInd,:),3)))) this formula is NOT for GRANDAVERAGE it's for matrix of ChanXtimeXepochs
        %double(squeeze(nanmean(sqrt(nansum([t(elecnum,timeInd,:)-nanmean(t(elecnum,timeInd,:),2)].^2,1)./size(elecnum,2)),3)));
        %    t=ALLEEG(i).data;
        % CH=[3 6 7 8 9]; LMGP=std(nanmean( squeeze(EEG.data(CH,:,:))   ,3)); % this is reza's LMFP - he say it's std of the channels
        LMFP(:,:,i)=smoothdata(LMFP(:,:,i),'gaussian',15);
        %double(mean(sqrt(sum((t(elecnum,timeInd,:)-mean(t(elecnum,timeInd,:),1)).^2)/size(t(elecnum,:,:),1)),3)); %Global Mean Field Power
        LMFPste(:,:,i)=double(squeeze(ste(std(t(elecnum,timeInd,:),0,1),3)));
        % double(squeeze(std(nanmean(t(elecnum,:,:),3)))); this formula is NOT for GRANDAVERAGE it's for matrix of ChanXtimeXepochs
        %double(squeeze(ste(sqrt(nansum([t(elecnum,timeInd,:)-nanmean(t(elecnum,timeInd,:),2)].^2,1)./size(elecnum,2)),3)));
        LMFPste(:,:,i)=smoothdata(LMFPste(:,:,i),'gaussian',15);
        %double(ste(sqrt(sum((t(elecnum,timeInd,:)-mean(t(elecnum,timeInd,:),1)).^2)/size(t(elecnum,:,:),1)),3));
        %     elseif strcmp(LMFPgrph,'on') & size(elecnum,2)==1
        %                 warning('LMFP can be computed only with 2 or more electrodes')
    end
    
    
    %TEP(:,:,i)=nanmean( t, 3);
    if strcmp(ISP,'on')
        [elecHOMOnum elecHOMOname]=homoElec(ALLEEG(i),elec);
        % [elecHOMOnum elecHOMOname]=homoElec(ALLEEG(1),{'fc4'});
        %ISPInd=ALLEEG(i).times>=statwin(1)+ISP_lag & ALLEEG(i).times<=statwin(2)+ISP_lag;
        ISPVecste(:,:,i)=double(ste(mean(t(elecHOMOnum,timeInd,:),1),3));
        ISPVecste(:,:,i)=smoothdata(ISPVecste(:,:,i),'gaussian',7);
        ISPVec1(:,:,i)=double(nanmean(mean(t(elecHOMOnum,timeInd,:),1),3));
        ISPVec1(:,:,i)=smoothdata(ISPVec1(:,:,i),'gaussian',7);
    end
    % trying to smooth stich after average
    %TEP(:,:,i)=nanmean(ALLEEG(i).data, 3);
    %TEP(:,stich,i)=hampel(TEP(:,stich,i)); %
    % TEP(:,stich,i)=sgolayfilt(double(TEP(:,stich,i)),3,11);
    %sgolayfilt(x,3,11)
    
    if strcmp(rectyfied,'on')
        TEPrecste(:,:,i)=double(ste(abs(mean(t(elecnum,timeInd,:),1)),3));
        TEPrecste(:,:,i)=smoothdata(TEPrecste(:,:,i),'gaussian',7);
        TEPrec(:,:,i)=double(nanmean(abs(mean(t(elecnum,timeInd,:),1)),3));
        TEPrec(:,:,i)=smoothdata(TEPrec(:,:,i),'gaussian',7);
        %TEPrec(:,:,i)=std(mean(t(elecnum,timeInd,:),1),0,3);
        
        ISPrecste(:,:,i)=double(ste(abs(mean(t(elecHOMOnum,timeInd,:),1)),3));
        ISPrecste(:,:,i)=smoothdata(ISPrecste(:,:,i),'gaussian',7);
        ISPrec(:,:,i)=double(nanmean(abs(mean(t(elecHOMOnum,timeInd,:),1)),3));
        ISPrec(:,:,i)=smoothdata(ISPrec(:,:,i),'gaussian',7);
        %ISPrec(:,:,i)=std(mean(t(elecHOMOnum,timeInd,:),1),0,3);
    end
    
    %TEPVecstat(:,:,i)=nanmean(TEP(:,statInd,i), 2);
    %figure;hold on; plot(EEG.times(timeInd),recVec1)
    %convbin1=[0 1 1 1 0];
    %ERPVec2=ERPVec1%conv(ERPVec1,convbin1,'same');
    
    clear statind timeind t stich
end
%export_fig diff_ADHD_Healthy_stp0.pdf -q101 -painters -append

%% wave PLOTS

if strcmp(GMFPgrph,'on')
    figure('position', [500 200 600 400]); hold on;
    f=fill([statwin statwin(2) statwin(1)],[amp amp 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
    patch([ALLEEG(dataset(:,1)).times(timeInd) fliplr(ALLEEG(dataset(:,1)).times(timeInd))],[GMFP(:,:,dataset(:,1))-GMFPste(:,:,dataset(:,1)) fliplr(GMFP(:,:,dataset(:,1))+GMFPste(:,:,dataset(:,1)))],[color{dataset(:,1)}.*0.7],'EdgeColor','none','FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %,'FaceVertexAlphaData',0.7,'FaceAlpha',0.3
    p1=plot(ALLEEG(dataset(:,1)).times(timeInd),GMFP(:,:,dataset(:,1)),'LineWidth',2);
    p1.Color=color{dataset(:,1)};
    patch([ALLEEG(dataset(:,2)).times(timeInd) fliplr(ALLEEG(dataset(:,2)).times(timeInd))],[GMFP(:,:,dataset(:,2))-GMFPste(:,:,dataset(:,2)) fliplr(GMFP(:,:,dataset(:,2))+GMFPste(:,:,dataset(:,2)))],[color{dataset(:,2)}.*0.7],'EdgeColor','none','FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.7,'FaceAlpha',0.3  ,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
    p2=plot(ALLEEG(dataset(:,2)).times(timeInd),GMFP(:,:,dataset(:,2)),'LineWidth',2);
    p2.Color=color{dataset(:,2)};
    axis([timeWin 0 amp]);
    %set(gca, 'ylimmode', 'auto')
    yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
    fb=fill([0 15 15 0],[yy yy 0 0],'w','edgecolor','none');
    ax=gca; ax.XAxisLocation='origin' ;  ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin';ax.XTick = [ -20 0 15 30 45 60 100 200 250];legend('show');clear title;% 0 15 30 60 80 100 120 140 160 180 200 250
    title(['Global Mean Field Power  _time:' strrep(num2str(statwin),'  ','_')]); ylabel ('\muv'); xlabel ('miliseconds');
    legend([p1 p2],{strrep([ALLEEG(dataset(:,1)).setname ' n=' num2str(size(ALLEEG(dataset(:,1)).subjects,2))],'_',' '),...
        strrep([ALLEEG(dataset(:,2)).setname ' n=' num2str(size(ALLEEG(dataset(:,2)).subjects,2))],'_',' ')})
    % ['ADHD  n=' num2str(size(corrsADHD(NnA,:),1))]
    
    hold off;
    clearvars f p1 p2 ax
    % stats
end

if strcmp(LMFPgrph,'on') %& size(elecnum,2)>1
    figfig1=figure('position', [500 200 600 400]); hold on;
    patch([ALLEEG(dataset(:,1)).times(timeInd) fliplr(ALLEEG(dataset(:,1)).times(timeInd))],[LMFP(:,:,dataset(:,1))-LMFPste(:,:,dataset(:,1)) fliplr(LMFP(:,:,dataset(:,1))+LMFPste(:,:,dataset(:,1)))],color{dataset(:,1)},'EdgeColor',[color{dataset(:,1)}.*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %
    p1=plot(ALLEEG(dataset(:,1)).times(timeInd),LMFP(:,:,dataset(:,1)),'LineWidth',2);
    p1.Color=color{dataset(:,1)};
    patch([ALLEEG(dataset(:,2)).times(timeInd) fliplr(ALLEEG(dataset(:,2)).times(timeInd))],[LMFP(:,:,dataset(:,2))-LMFPste(:,:,dataset(:,2)) fliplr(LMFP(:,:,dataset(:,2))+LMFPste(:,:,dataset(:,2)))],color{dataset(:,2)},'EdgeColor',[color{dataset(:,2)}.*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
    p2=plot(ALLEEG(dataset(:,2)).times(timeInd),LMFP(:,:,dataset(:,2)),'LineWidth',2);
    p2.Color=color{dataset(:,2)};
    yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
    f=fill([statwin statwin(2) statwin(1)],[yy yy 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
    chH = get(gca,'Children');
    set(gca,'Children',[chH(2:end) ;chH(1)])
    fb=fill([0 15 15 0],[yy yy 0 0],'w','edgecolor','none');
    %axis([timeWin 0 amp]);
    ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
    legend('show'); clear title; title(['Local Mean Field Power ' strjoin([elecname]) '  time:' strrep(num2str(statwin),'  ','_')], 'interpreter', 'none'); ylabel ('\muv'); xlabel ('miliseconds');
    legend([p1 p2],{strrep([ALLEEG(dataset(:,1)).setname ' n=' num2str(size(ALLEEG(dataset(:,1)).subjects,2))],'_',' '),...
    strrep([ALLEEG(dataset(:,2)).setname ' n=' num2str(size(ALLEEG(dataset(:,2)).subjects,2))],'_',' ')})
    set(ax ,'Layer', 'Top')
    % ['ADHD  n=' num2str(size(corrsADHD(NnA,:),1))]
    hold off;
    clearvars f p1 p2 ax
    % stats
    % elseif strcmp(LMFPgrph,'on') & size(elecnum,2)==1
    %                 warning('LMFP can be computed only with 2 or more electrodes')
end

if strcmp(wavgrph,'on')
    
    figure('position', [500 200 600 400]); hold on;
    f=fill([statwin statwin(2) statwin(1)],[amp-0.5 amp-0.5 -(amp-0.5) -(amp-0.5)],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
    timeInd=ALLEEG(1).times>=timeWin(1) & ALLEEG(1).times<=timeWin(2);
    patch([ALLEEG(dataset(:,1)).times(timeInd) fliplr(ALLEEG(dataset(:,1)).times(timeInd))],...
        [TEPVec1(:,:,dataset(:,1))-TEPVecste(:,:,dataset(:,1)) fliplr(TEPVec1(:,:,dataset(:,1))+TEPVecste(:,:,dataset(:,1)))],color{dataset(:,1)},'FaceVertexAlphaData',0.1,'LineWidth',0.2,'FaceAlpha',0.3,'EdgeColor',[color{dataset(:,1)}.*0.5]) %,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3,'EdgeColor',color{dataset(:,1)}
    p1=plot(ALLEEG(dataset(:,1)).times(timeInd),TEPVec1(:,:,dataset(:,1)),'LineWidth',2.5);
    p1.Color=color{dataset(:,1)};
    timeInd=ALLEEG(2).times>=timeWin(1) & ALLEEG(2).times<=timeWin(2);
    patch([ALLEEG(dataset(:,2)).times(timeInd) fliplr(ALLEEG(dataset(:,2)).times(timeInd))],...
        [TEPVec1(:,:,dataset(:,2))-TEPVecste(:,:,dataset(:,2)) fliplr(TEPVec1(:,:,dataset(:,2))+TEPVecste(:,:,dataset(:,2)))],color{dataset(:,2)},'FaceVertexAlphaData',0.1,'LineWidth',0.2,'FaceAlpha',0.3,'EdgeColor',[color{dataset(:,2)}.*0.5])%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3,'EdgeColor',color{dataset(:,2)}
    p2=plot(ALLEEG(dataset(:,2)).times(timeInd),TEPVec1(:,:,dataset(:,2)),'LineWidth',2.5);
    p2.Color=color{dataset(:,2)};
    yy=max(get(gca,'YLim')); 
    yymin=min(get(gca,'YLim'));
    fb=fill([0 25 25 0],[yy yy yymin yymin],'w','edgecolor','none');
    axis([timeWin -amp amp]);
    ax=gca;  ax.FontName='Helvetica Neue' ; ax.FontSize=11;  ax.XAxisLocation='origin' ;
    ax.YAxisLocation = 'origin';    ax.XTick = [ -20 0 15 30 45 60 100 200 250];
    ax.FontSize=11; legend('show');clear title; %ax.XTick = [timeWin(1):100:timeWin(2)]
    title([strjoin([elecname]) strrep(num2str(statwin),'  ','_')], 'interpreter', 'none'); ylabel ('\muv'); xlabel ('miliseconds');
    legend([p1 p2],{strrep([ALLEEG(dataset(:,1)).setname ' n=' num2str(size(ALLEEG(dataset(:,1)).subjects,2))],'_',' '),...
        strrep([ALLEEG(dataset(:,2)).setname ' n=' num2str(size(ALLEEG(dataset(:,2)).subjects,2)) ],'_',' ')})
    hold off;
    clearvars f p1 p2 ax
    % stats
end

if strcmp(ISP,'on')
    figure('position', [500 200 600 400]); hold on;
    f=fill([statwin statwin(2) statwin(1)]+ISP_lag,[amp-0.5 amp-0.5 -(amp-0.5) -(amp-0.5)],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
    p1=plot(ALLEEG(dataset(:,1)).times(timeInd),ISPVec1(:,:,dataset(:,1)),'LineWidth',2);
    p1.Color=color{dataset(:,1)};
    patch([ALLEEG(dataset(:,1)).times(timeInd) fliplr(ALLEEG(dataset(:,1)).times(timeInd))],[ISPVec1(:,:,dataset(:,1))-ISPVecste(:,:,dataset(:,1)) fliplr(ISPVec1(:,:,dataset(:,1))+ISPVecste(:,:,dataset(:,1)))],color{dataset(:,1)},'EdgeColor',[color{dataset(:,1)}.*0.5],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
    p2=plot(ALLEEG(dataset(:,2)).times(timeInd),ISPVec1(:,:,dataset(:,2)),'LineWidth',2);
    p2.Color=color{dataset(:,2)};
    patch([ALLEEG(dataset(:,2)).times(timeInd) fliplr(ALLEEG(dataset(:,2)).times(timeInd))],[ISPVec1(:,:,dataset(:,2))-ISPVecste(:,:,dataset(:,2)) fliplr(ISPVec1(:,:,dataset(:,2))+ISPVecste(:,:,dataset(:,2)))],color{dataset(:,2)},'EdgeColor',[color{dataset(:,2)}.*0.5],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
    axis([timeWin -amp amp]);
    ax=gca; ax.XAxisLocation='origin' ;  ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin';ax.XTick = [-20 0 15 30 45 60 100 200 250];legend('show');clear title;
    title(['Homologic ISP ' strjoin([elecHOMOname])], 'interpreter', 'none'); ylabel ('\muv'); xlabel ('miliseconds');
    legend([p1 p2],{strrep([ALLEEG(dataset(:,1)).setname(1:end) ' n=' num2str(size(ALLEEG(dataset(:,1)).subjects,2))],'_',' '),...
        strrep([ALLEEG(dataset(:,2)).setname(1:end) ' n=' num2str(size(ALLEEG(dataset(:,2)).subjects,2)) ],'_',' ')})
    hold off;
    clearvars f p1 p2 ax
end

if strcmp(rectyfied,'on')
    
    figure('position', [500 200 600 400]); hold on;
    f=fill([statwin statwin(2) statwin(1)],[amp-0.5 amp-0.5 -(amp-0.5) -(amp-0.5)],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
    p1=plot(ALLEEG(dataset(:,1)).times(timeInd),TEPrec(:,:,dataset(:,1)),'LineWidth',2);
    p1.Color=color{dataset(:,1)};
    patch([ALLEEG(dataset(:,1)).times(timeInd) fliplr(ALLEEG(dataset(:,1)).times(timeInd))],...
        [TEPrec(:,:,dataset(:,1))-TEPrecste(:,:,dataset(:,1)) fliplr(TEPrec(:,:,dataset(:,1))+TEPrecste(:,:,dataset(:,1)))],color{dataset(:,1)},'EdgeColor',[color{dataset(:,1)}.*0.5],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)% ,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
    p2=plot(ALLEEG(dataset(:,2)).times(timeInd),TEPrec(:,:,dataset(:,2)),'LineWidth',2);
    p2.Color=color{dataset(:,2)};
    patch([ALLEEG(dataset(:,2)).times(timeInd) fliplr(ALLEEG(dataset(:,2)).times(timeInd))],...
        [TEPrec(:,:,dataset(:,2))-TEPrecste(:,:,dataset(:,2)) fliplr(TEPrec(:,:,dataset(:,2))+TEPrecste(:,:,dataset(:,2)))],color{dataset(:,2)},'EdgeColor',[color{dataset(:,2)}.*0.5],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
    axis([timeWin -amp amp]);
    ax=gca; ax.XAxisLocation='origin' ;  ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin';ax.XTick = [-20 0 15 30 45 60 100 200 250];legend('show');clear title;
    title(['Rectyfied ' strjoin([elecname])], 'interpreter', 'none'); ylabel ('\muv'); xlabel ('miliseconds');
    legend([p1 p2],{strrep([ALLEEG(dataset(:,1)).setname ' n=' num2str(size(ALLEEG(dataset(:,1)).subjects,2))],'_',' '),...
        strrep([ALLEEG(dataset(:,2)).setname ' n=' num2str(size(ALLEEG(dataset(:,2)).subjects,2)) ],'_',' ')})
    hold off;
    clearvars f p1 p2 ax
    
    if strcmp(ISP,'on')
        figure('position', [500 200 600 400]); hold on;
        f=fill([statwin statwin(2) statwin(1)]+ISP_lag,[amp-0.5 amp-0.5 -(amp-0.5) -(amp-0.5)],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
        p1=plot(ALLEEG(dataset(:,1)).times(timeInd),ISPrec(:,:,dataset(:,1)),'LineWidth',2);
        p1.Color=color{dataset(:,1)};
        patch([ALLEEG(dataset(:,1)).times(timeInd) fliplr(ALLEEG(dataset(:,1)).times(timeInd))],[ISPrec(:,:,dataset(:,1))-ISPrecste(:,:,dataset(:,1)) fliplr(ISPrec(:,:,dataset(:,1))+ISPrecste(:,:,dataset(:,1)))],color{dataset(:,1)},'EdgeColor',[color{dataset(:,1)}.*0.5],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
        p2=plot(ALLEEG(dataset(:,2)).times(timeInd),ISPrec(:,:,dataset(:,2)),'LineWidth',2);
        p2.Color=color{dataset(:,2)};
        patch([ALLEEG(dataset(:,2)).times(timeInd) fliplr(ALLEEG(dataset(:,2)).times(timeInd))],[ISPrec(:,:,dataset(:,2))-ISPrecste(:,:,dataset(:,2)) fliplr(ISPrec(:,:,dataset(:,2))+ISPrecste(:,:,dataset(:,2)))],color{dataset(:,2)},'EdgeColor',[color{dataset(:,2)}.*0.5],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
        axis([timeWin -amp amp]);
        ax=gca; ax.XAxisLocation='origin' ;  ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin';ax.XTick = [ 0 15 30 60 100 200 250];
        ax.FontSize=11; legend('show');clear title;
        title(['Homologic rectifyed ISP ' strjoin([elecHOMOname])], 'interpreter', 'none'); ylabel ('\muv'); xlabel ('miliseconds');
        legend([p1 p2],{strrep([ALLEEG(dataset(:,1)).setname ' n=' num2str(size(ALLEEG(dataset(:,1)).subjects,2))],'_',' '),...
            strrep([ALLEEG(dataset(:,2)).setname ' n=' num2str(size(ALLEEG(dataset(:,2)).subjects,2)) ],'_',' ')})
        hold off;
        clearvars f p1 p2 ax
        %subs.ISP=(permute(sum(abs(mean(ALLEEG(i).data(elecHOMOnum,ISPInd,:),1)),2), [3 2 1])./permute(sum(abs(mean(ALLEEG(i).data(elecnum,statInd,:),1)),2), [3 2 1])).*100;
    end
end
%% stats

warning('off')

%dat= [3 4 5 6]
% figure;  hold on;
for i=dataset %size({ALLEEG.setname},2)
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, find(strcmp({ALLEEG.setname},EEG.setname)),'retrieve',i,'study',0);
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',CURRENTSTUDY);
    %elecnum=elecName(EEG,elec);
    %[elecnum elecname]=elecName(ALLEEG(i),elec);
    subsT=table;
    subsT.Subjects=([ALLEEG(i).subjects]'); % str2double if the subjects field is string
    %subs.Subjects=[ALLEEG(i).subjects]' % if the subjects are double
    Group=ALLEEG(i).filename;%(1:end);
    if ~isempty(strfind(lower(Group),lower(Groups{1})))
        Group=Groups{1};
    elseif ~isempty(strfind(lower(Group),lower(Groups{2})))
        Group=Groups{2};
    end
    %cond=ALLEEG(i).filename(strfind(ALLEEG(i).filename,condName):strfind(ALLEEG(i).filename,condName)+3);
    subsT.group=repmat({Group},size(subsT,1),1);
    subsT.group=categorical(subsT.group);
    statInd=ALLEEG(i).times>=statwin(1) & ALLEEG(i).times<=statwin(2);
    timeInd=ALLEEG(i).times>=timeWin(1) & ALLEEG(i).times<=timeWin(2);
    
    %head averages
    for s=1:size(subsT,1)
        head{i}(:,s)=double(mean(ALLEEG(i).data(:,timeInd,s), 2));
        % head{i}=heads;
    end
    
    switch method
        case 'mean' %strcmp(method,'mean')
            subsT.stat=double(squeeze(mean(mean(ALLEEG(i).data(elecnum,statInd,:),1),2)));
            %head=[];
            % head{i}(:,:,s)=mean(ALLEEG(i).data(:,statInd,s), 2);
        case 'std' %strcmp(method,'std')
            subsT.stat=double(squeeze(std(mean(ALLEEG(i).data(elecnum,statInd,:),1),1,2)));
            %head=[];
        case 'rectified' %strcmp(method,'rectified')
            subsT.stat=double(squeeze(abs(mean(mean(ALLEEG(i).data(elecnum,statInd,:),1),2))));
            %head=[];
        case 'AUC'
            subsT.stat=double(squeeze(sum(abs(mean(ALLEEG(i).data(elecnum,statInd,:),1)),2)))./(ALLEEG(i).srate/1000);
            %Int = trapz(X, Y)
            %for sub2=1:size(ALLEEG(i).data,3)
            %subsT.stat(sub)=double(trapz(ALLEEG(i).times(statInd),abs(mean(ALLEEG(i).data(elecnum,statInd,sub),1))));
            %subsT.stat(sub2)=double(sum(abs(mean(ALLEEG(i).data(elecnum,statInd,sub2),1))))./(ALLEEG(i).srate/1000);
            %end
            %clear sub
            %head=[];
        case 'GMFP'
            %elecnum_G=size(ALLEEG(i).data,1);
            subsT.stat=double(squeeze(nanmean(std(ALLEEG(i).data(:,statInd,:),0,1),2)));
            %double(squeeze(nanmean(sqrt(nansum([ALLEEG(i).data([1:elecnum_G],statInd,:)-nanmean(ALLEEG(i).data([1:elecnum_G],statInd,:),2)].^2,1)./size([1:elecnum_G],2)),2)));
            %double(squeeze(mean(sqrt(sum((ALLEEG(i).data(:,statInd,:)-mean(ALLEEG(i).data(:,statInd,:),1)).^2)/size(ALLEEG(i).data,1)),2)));
            %head=[];
            % gmpf(i,j,k,ALLEEG(ind).subjects(p)-100,q) =sqrt(sum((temp(:,q,p)-mean(temp(:,q,p),1)).^2)/63);%Global Mean Field Power - sam's quation
        case 'LMFP'
            subsT.stat=double(squeeze(nanmean(std(ALLEEG(i).data(elecnum,statInd,:),0,1),2)));
            %std(nanmean(ALLEEG(i).data(elecnum,:,:),3));
            %double(squeeze(nanmean(sqrt(nansum([ALLEEG(i).data(elecnum,statInd,:)-nanmean(ALLEEG(i).data(elecnum,statInd,:),2)].^2,1)./size(elecnum,2)),2)));
            %   double(squeeze(nanmean(sqrt(nansum([ALLEEG(i).data(elecnum,:,:)-nanmean(ALLEEG(i).data(elecnum,:,:),1)].^2,1)./size(elecnum,2)),2)));
            %             double(squeeze(mean(sqrt(sum((ALLEEG(i).data(:,statInd,:)-mean(ALLEEG(i).data(:,statInd,:),1)).^2)/size(ALLEEG(i).data,1)),2)));   timeInd
            %             if size(elecnum,2)>1
            %             subsT.stat=double(squeeze(mean(sqrt(sum((ALLEEG(i).data(elecnum,statInd,:)-mean(ALLEEG(i).data(elecnum,statInd,:),1)).^2)/size(ALLEEG(i).data(elecnum,:,:),1)),2)));
            %             else
            %                 error('LMFP can be computed only with 2 or more electrodes')
            %             end
        case 'admean' %strcmp(method,'admean')%elseif strcmp(method,'admean')
            
            for s=1:size(subsT,1)
                meanvec=double(mean(ALLEEG(i).data(elecnum,:,s), 1));
                % tmplot(s,:)=meanvec(statInd);
                nanvec=NaN(1, size(meanvec,2));
                nanvec(statInd)=meanvec(statInd);
                if strcmp(peaktype, 'positive')
                    [vmax, ind]=max(nanvec);
                elseif strcmp(peaktype, 'negative')
                    [vmax, ind]=min(nanvec);
                end
                timeind=ALLEEG(i).times>=ALLEEG(i).times(ind)-smallWinHalfLength & ALLEEG(i).times<=ALLEEG(i).times(ind)+smallWinHalfLength;
                %head=[];
                % head{i}(:,s)=mean(ALLEEG(i).data(:,timeind,s), 2);
                %headminmax{i}(:,:,s)=max(ALLEEG(i).data(:,timeInd,s),[],2)-min(ALLEEG(i).data(:,timeInd,s),[],2);
                if iscell(subsT.Subjects) & ~iscell(ALLEEG(i).subjects(:)) %strcmp(class(subsT.Subjects),'cell')
                    subsT.Subjects=str2num(cell2mat(subsT.Subjects)); % str2num(cell2mat(subsT.Subjects));
                end
                [~, tt, ~]=intersect(subsT.Subjects,ALLEEG(i).subjects(s));
                subsT.stat(tt)=double(mean(meanvec(timeind)));
                
                
                clearvars timeind meanvec ind vmax nanvec
            end
            
        otherwise
            warning ('optional arrguments: mean / std / rectified / AUC / GMFP / LMFP / admean')
            %   plot(ALLEEG(i).times(statInd),tmplot) %plotting all subjects TEP segment
    end
    
    if iscell(subsT.Subjects) %strcmp(class(subsT.Subjects),'cell')
        subsT.Subjects=cell2mat(subsT.Subjects); % str2num(cell2mat(subsT.Subjects));
    end
    if isstring(subsT.Subjects)
        subsT.Subjects=str2num(subsT.Subjects);
    end
    
    
    if strcmp(method,'admean')
        subsT.Properties.VariableNames{'stat'} = strrep(strjoin({name method strrep(strjoin(elecname),' ','_') strrep(num2str([statwin(1) smallWinHalfLength statwin(2)]),'  ','_')}, '_'),'__','_');
    elseif strcmp(method,'GMFP')
        subsT.Properties.VariableNames{'stat'} = strjoin({name method strrep(num2str(statwin),'  ','_')}, '_');
    else
        subsT.Properties.VariableNames{'stat'} = strjoin({name method strrep(strjoin(elecname),' ','_') strrep(num2str(statwin),'  ','_')}, '_');
    end
    
    if strcmp(ISP,'on')
        ISPInd=ALLEEG(i).times>=statwin(1)+ISP_lag & ALLEEG(i).times<=statwin(2)+ISP_lag;
        
        %subsT.ISP=(permute(sum(abs(mean(ALLEEG(i).data(elecHOMOnum,ISPInd,:),1)),2), [3 2 1])./permute(sum(abs(mean(ALLEEG(i).data(elecnum,statInd,:),1)),2), [3 2 1])).*100;
        for sub=1:size(ALLEEG(i).data,3)
            %subsT.stat(sub)=trapz(ALLEEG(i).times(statInd),abs(mean(ALLEEG(i).data(elecnum,statInd,sub),1)));
            subsT.ISP(sub)= double(trapz(ALLEEG(i).times(ISPInd),abs(mean(ALLEEG(i).data(elecHOMOnum,ISPInd,sub),1)))./trapz(ALLEEG(i).times(statInd),abs(mean(ALLEEG(i).data(elecnum,statInd,sub),1))));
        end
        subsT.Properties.VariableNames{'ISP'} = strrep(strjoin({'ISP' strrep(strjoin([elecname elecHOMOname]),' ','_') strrep(num2str([statwin(1) statwin(2)]),'  ','_')}, '_'),'__','_');
        
    end
    if iscell(subsT.Subjects)%(class(subsT.Subjects),'cell')
        subsT.Subjects=str2num(cell2mat(subsT.Subjects)); % str2num(cell2mat(subsT.Subjects));
    end
    subsT=sortrows(subsT,1);
    stra{i}=table2struct(subsT);
    %writetable(subsT,[name '_' cond '_' Group ' timewin-' num2str(statwin) '_' strjoin(elec) '.csv']);
    
end


for i=1:size(stra,2)
    if isempty(stra{i})
    else
        f(i,:)=fieldnames(stra{i});
        fi(:,i)=f(i,3);
    end
end
fi=unique(fi(~cellfun(@isempty, fi)));

for i=1:size(stra,2)
    if isempty(stra{i})
    else
        ind(i)=double(nansum(strcmp(fieldnames(stra{i}),fi(1))));
    end
end

x=find(ind);
a=struct2table(vertcat(stra{1,x(1)},stra{1,x(2)}));
%a.group=categorical(a.group); %line added at 20170222
stats=a;
%stat

if strcmpi(rem_outlie,'on')
    clear remout1 remout2
    remout1=stats{stats.group==Groups(1),3};
    remout1(isoutlier(remout1))=nan;
    stats{stats.group==Groups(1),3}=remout1;
    remout2=stats{stats.group==Groups(2),3};
    remout2(isoutlier(remout2))=nan;
    stats{stats.group==Groups(2),3}=remout2;
    
end

% if iscell(stats.Subjects)%(class(subsT.Subjects),'cell')
%     stats.Subjects=str2num(cell2mat(stats.Subjects)); % str2num(cell2mat(subsT.Subjects));
% elseif ischar(stats.Subjects)
%     stats.Subjects=str2num(stats.Subjects);
% end

stats

[~, P , ~, desc]=ttest2(double(stats.(stats.Properties.VariableNames{3})(stats.group==Groups(2))),double(stats.(stats.Properties.VariableNames{3})(stats.group==Groups(1))))
%P=signrank(double(stats.(stats.Properties.VariableNames{3})(stats.group==Groups(2))),double(stats.(stats.Properties.VariableNames{3})(stats.group==Groups(1))))
%P=ranksum(double(stats.(stats.Properties.VariableNames{3})(stats.group==Groups(2))),double(stats.(stats.Properties.VariableNames{3})(stats.group==Groups(1))))

if strcmp(ISP,'on')
    % [~, P_isp]=ttest2(double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(2))),double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(1))))
    %P_isp=signrank(double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(2))),double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(1))))
    %P_isp=ranksum(double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(2))),double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(1))))
    [~, P_isp , ~, isp_desc]=ttest2(double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(2))),double(stats.(stats.Properties.VariableNames{4})(stats.group==Groups(1))))
end
% repeated ANOVA - not necessery here
% model=[stats.Properties.VariableNames{3} '-' stats.Properties.VariableNames{4} '~' stats.Properties.VariableNames{2}];
% repeated=[3 5];
% rm = fitrm(stats, model, 'WithinDesign', repeated);
% no = ranova(rm)
% post = multcompare(rm,'group')
head=  head( find(~cellfun(@isempty,head)));


%% stats Group plot
if strcmp(bargrph,'on')
    
    figure('position', [200 300 300 300]);hold on;
    
    bar(1,nanmean(stats.(3)(stats.group==Groups(2))),'FaceColor',Healthycolor)
    bar(2,nanmean(stats.(3)(stats.group==Groups(1))),'FaceColor',ADHDcolor)
    set(gca,'XTick',[1 2],'XTickLabel',[[Groups(2)]  [Groups(1)]],'FontSize',13)
    %axis([0 0.45 -4 0])%-100*10^3 100*10^3])
    %annotation('textbox','String',['p=' num2str(p)]);
    errorbar([nanmean(stats.(3)(stats.group==Groups(2))) nanmean(stats.(3)(stats.group==Groups(1)))],...
        [ste(stats.(3)(stats.group==Groups(2))) ste(stats.(3)(stats.group==Groups(1)))], '+k' )
    %axis([0 3 -amp 0 ])
    legend(strrep(strjoin([Groups(2) ' n=' num2str(size(stats(stats.group==Groups(2),:),1))]),'_',' '),...
        strrep(strjoin([Groups(1) ' n=' num2str(size(stats(stats.group==Groups(1),:),1))]),'_',' '))
    
    title(['p=' num2str(P) newline strrep(strjoin(stats.Properties.VariableNames(3)),'_',' ')]); ylabel ('Amplitude (\muv)');
    hold off;
    
    if strcmp(ISP,'on')
        figure('position', [200 300 300 300]);hold on;
        
        bar(1,nanmean(stats.(4)(stats.group==Groups(2))),'FaceColor',Healthycolor)
        bar(2,nanmean(stats.(4)(stats.group==Groups(1))),'FaceColor',ADHDcolor)
        set(gca,'XTick',[1 2],'XTickLabel',[[Groups(2)]  [Groups(1)]],'FontSize',13)
        %axis([0 0.45 -4 0])%-100*10^3 100*10^3])
        %annotation('textbox','String',['p=' num2str(p)]);
        errorbar([nanmean(stats.(4)(stats.group==Groups(2))) nanmean(stats.(4)(stats.group==Groups(1)))],...
            [ste(stats.(4)(stats.group==Groups(2))) ste(stats.(4)(stats.group==Groups(1)))], '+k' )
        legend(strrep(strjoin([Groups(2) ' n=' num2str(size(stats(stats.group==Groups(2),:),1))]),'_',' '),...
            strrep(strjoin([Groups(1) ' n=' num2str(size(stats(stats.group==Groups(1),:),1))]),'_',' ') )
        title(['p=' num2str(P_isp) newline strrep(strjoin(stats.Properties.VariableNames(4)),'_',' ')]); ylabel ('AUC');
        hold off;
        
    end
    
end

if strcmp(violinplt,'on')
    
    ya=stats{stats.group==Groups(1),3};
    ya=ya(~isnan(ya));
    yb=stats{stats.group==Groups(2),3};
    yb=yb(~isnan(yb));
    
    colorsg=[ADHDcolor ; Healthycolor] ;
    
    y=[ya;yb];
    
    x=[cellstr(repmat([strrep(Groups{1},'_',' ') '-' strrep(Groups{2},'_',' ')],size(y,1),1))];
    group=[cellstr(repmat(strrep(Groups{1},'_',' '),size(ya,1),1)) ; cellstr(repmat(strrep(Groups{2},'_',' '),size(yb,1),1))];
    clear g
    g=gramm('x',x,'y',y,'color',group,'group',group);%,'color',bo,
    g.stat_violin('fill','transparent' ,'width',1 ,'dodge',0.5 );%'normalization','count', 'npoints',15,
    %g.stat_boxplot('width', 0.6 )
    g.stat_summary('type','sem' ,'geom','lines','dodge', 0.5 );%'quartile' '95percentile' 'std' 'fitnormalci'
    %g.geom_jitter('width', 0.2,'dodge', 0.5 )
    %g.geom_vline('xintercept',1)
    g.geom_polygon('line_style','k--')
    g.set_color_options('map',colorsg,'n_color',2,'n_lightness',1);
    g.set_text_options('font','Helvetica Neue','interpreter', 'tex',...
        'legend_title_scaling', 1.1 , 'big_title_scaling', 1.2 );
    g.set_names('x','','y',['Amplitude (' '\mu' 'v)']);
    g.set_title(['p=' num2str(P) '      t_(' '_{' num2str(desc.df) '}_)=' num2str(desc.tstat) newline strrep(strjoin(stats.Properties.VariableNames(3)),'_',' ')]);
    figure('Position',[100 100 320 450]);
    %g.coord_flip();
    %g.set_order_options('color',1);
    g.draw();
    %aaa=gca aaa.Children LineWidth=0.5; bbb=findall(aaa.Children, 'Type', 'Line');
    %bbb(4).LineWidth=3.2
    % g.export('file_name','MST_violin_pre-post_Hipp_SCD','file_type','pdf')
    
    if strcmp(ISP,'on')
        isp_ya=stats{stats.group==Groups(1),4};
        isp_ya=isp_ya(~isnan(isp_ya));
        isp_yb=stats{stats.group==Groups(2),4};
        isp_yb=isp_yb(~isnan(isp_yb));
        
        colorsg=[ADHDcolor ; Healthycolor] ;
        
        isp_y=[isp_ya;isp_yb];
        
        x=[cellstr(repmat([strrep(Groups{1},'_',' ') '-' strrep(Groups{2},'_',' ')],size(isp_y,1),1))];
        group=[cellstr(repmat(strrep(Groups{1},'_',' '),size(isp_ya,1),1)) ; cellstr(repmat(strrep(Groups{2},'_',' '),size(isp_yb,1),1))];
        clear g
        isp_g=gramm('x',x,'y',isp_y,'color',group,'group',group);%,'color',bo,
        isp_g.stat_violin('fill','transparent' ,'width',1 ,'dodge',0.5 );%'normalization','count', 'npoints',15,
        %g.stat_boxplot('width', 0.6 )
        isp_g.stat_summary('type','sem' ,'geom','lines','dodge', 0.5 );%'quartile' '95percentile' 'std' 'fitnormalci'
        %g.geom_jitter('width', 0.2,'dodge', 0.5 )
        %g.geom_vline('xintercept',1)
        isp_g.set_color_options('map',colorsg,'n_color',2,'n_lightness',1);
        isp_g.set_text_options('font','Helvetica Neue','interpreter', 'tex',...
            'legend_title_scaling', 1.1 , 'big_title_scaling', 1.2 );
        isp_g.set_names('x','','y',['Amplitude (' '\mu' 'v)']);
        isp_g.set_title(['p=' num2str(P_isp) '      t_(' '_{' num2str(isp_desc.df) '}_)=' num2str(isp_desc.tstat) newline strrep(strjoin(stats.Properties.VariableNames(4)),'_',' ')]);
        figure('Position',[100 100 320 450]);
        %g.coord_flip();
        %g.set_order_options('color',1);
        isp_g.draw();
        
    end
end
warning('on')
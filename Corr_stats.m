function [corMatHealthyF,corMatADHDF, corMatF, ff] = Corr_stats(Final, corrvars, cntrgroup, expgroup, alpha, cortype, disregard_zeros, disregard_num ,rem_outlie) 
% corrvars=[3 4] corrvars=[105 107] corrvars=[55 57]  corrvars=[18 30] 
%Final2=Final corrvars=[21 64];corrvars=[141 153];
% Final=Final2 % Final=Final_pre_post; % Final=Final3
% cortype='pearson'
%cntrgroup= 'Healthy' ; expgroup='MDD'. expgroup='ADHD'.
%cntrgroup= 'sham' ; expgroup='real'
% alpha=0.05
% disregard_zeros='on'
% rem_outlie='on'

% corrvars=[8 20];
% alpha=0.05;
% disregard_zeros='on';
% rem_outlie='on';
% cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
% 
% [corMatADHDF, corMatHealthyF, corMatF] = Corr_stats(Final2,vars, 'Healthy', 'MDD', alpha, cortype,disregard_zeros, rem_outlie );
% 
%expgroup='real';  cntrgroup='sham';
warning('off')
if sum(strcmpi(Final.Properties.VariableNames,'filt'))>0
    Final=Final(Final.filt==1,:);
end


vars=corrvars;
groups=Final.Group;
% corrsADHD=table2array(Final(Final.group==expgroup,3:end)) ;
tabdbl=~strcmp(varfun(@class,Final,'OutputFormat','cell'),'double') & ~strcmp(varfun(@class,Final,'OutputFormat','cell'),'single');

fin2=table2cell(Final); fin2(:,tabdbl)=num2cell(nan(size(Final,1), sum(tabdbl)));
fin2(:,strcmp(varfun(@class,Final,'OutputFormat','cell'),'single'))=num2cell(double(cell2mat(fin2(:,strcmp(varfun(@class,Final,'OutputFormat','cell'),'single')))));
corrsADHD=cell2mat(fin2(groups==expgroup,:));
%corrsADHD=table2array(Final(Final.group==expgroup,:)) ;
corrsHealthy=cell2mat(fin2(groups==cntrgroup,:));
%corrsHealthy=table2array(Final(Final.group==cntrgroup,3:end));
corrs=cell2mat(fin2);
%corrs=table2array(Final(:,3:end));

% zeros removal
if strcmpi(disregard_zeros,'on')
    corrsADHD(corrsADHD==0)=nan;
    corrsHealthy(corrsHealthy==0)=nan;
    corrs(corrs==0)=nan;
end

if disregard_num~=0
    corrsADHD(corrsADHD==disregard_num)=nan;
    corrsHealthy(corrsHealthy==disregard_num)=nan;
    corrs(corrs==disregard_num)=nan;
end

%outlier removal
if strcmpi(rem_outlie,'on')
    corrsADHD(isoutlier(corrsADHD))=nan;
    corrsHealthy(isoutlier(corrsHealthy))=nan;
    corrs(isoutlier(corrs))=nan;
end

%tab_var_nam=strcat('Var_', sprintfc('%d',[3:size(Final,2)-2]),'_', Final.Properties.VariableNames(3:end));
tab_var_nam=Final.Properties.VariableNames;
% ADHD
%[cormatADHD.R cormatADHD.P]=corrcoef(corrsADHD,'rows','pairwise'); %'rows','pairwise'
%[cormatADHD.R, cormatADHD.P] = corr(Final{groups==expgroup,vars},'rows','pairwise', 'type', cortype);%'rows','pairwise'
[cormatADHD.R, cormatADHD.P] = corr(corrsADHD,'rows','pairwise', 'type', cortype);%'rows','pairwise'


corMatADHD=nan(size(cormatADHD.R));
corMatADHD(cormatADHD.P<alpha)=cormatADHD.R(cormatADHD.P<alpha);
corMatADHDF=array2table(corMatADHD,'VariableNames',tab_var_nam,'RowNames',tab_var_nam');
%corMatADHDF=corMatADHDF(sum(~isnan(corMatADHDF{:,:}),2)>=1,sum(~isnan(corMatADHDF{:,:}))>=1);


% Healthy
%[cormatHealthy.R cormatHealthy.P]=corrcoef(corrsHealthy,'rows','pairwise');
%[cormatHealthy.R, cormatHealthy.P]=corr(Final{groups==cntrgroup,vars},'rows','pairwise','type', cortype); %'rows','pairwise'
[cormatHealthy.R, cormatHealthy.P]=corr(corrsHealthy,'rows','pairwise','type', cortype); %'rows','pairwise'

corMatHealthy=nan(size(cormatHealthy.R));
corMatHealthy(cormatHealthy.P<alpha)=cormatHealthy.R(cormatHealthy.P<alpha);
corMatHealthyF=array2table(corMatHealthy,'VariableNames',tab_var_nam,'RowNames',tab_var_nam');
%corMatHealthyF=corMatHealthyF(sum(~isnan(corMatHealthyF{:,:}),2)>=1,sum(~isnan(corMatHealthyF{:,:}))>=1);

% general
%[cormat.R, cormat.P]=corr(Final{:,vars},'rows','pairwise', 'type', cortype);%'rows','pairwise'
%[cormat.R cormat.P]=corrcoef(corrs,'rows','pairwise');
[cormat.R cormat.P]=corr(corrs,'rows','pairwise','type', cortype);
%corrplot(corrs,'rows','pairwise', 'type', cortype)
corMat=nan(size(cormat.R));
corMat(cormat.P<alpha)=cormat.R(cormat.P<alpha);
%corMat=cormat.R.*(cormat.P<alpha);
corMatF=array2table(corMat,'VariableNames',tab_var_nam,'RowNames',tab_var_nam');
%corMatF=corMatF(sum(~isnan(corMatF{:,:}),2)>=1,sum(~isnan(corMatF{:,:}))>=1);
% corMatF2=corMatF
%% PLOT VARS correlation
co=brewermap(4,'PuOr'); %PuOr RdBu
ADHDcolor=co(2,:);
ADHDcolor1=co(1,:);%'r'; %[0 0 0];
Healthycolor=co(3,:); % 'b' ;%[166/255 166/255 166/255];
Healthycolor1=co(4,:);
%ADHDcolor=[238/255 17/255 17/255];%'r'; %[0 0 0];
%Healthycolor=[43/255 87/255 151/255]; % 'b' ;%[166/255 166/255 166/255];
%ADHDcolor=[0 0 0]; 
%Healthycolor=[166/255 166/255 166/255]; 
%vars=[16 68]
%Nn= ~isnan(corrs(:,vars(1))) & ~isnan(corrs(:,vars(2)));
Nn= isfinite(corrs(:,vars(1))) & isfinite(corrs(:,vars(2)));
%cortype='Spearman';
%[r p]=corr(corrs(Nn,vars(1)),corrs(Nn,vars(2)),'rows', 'complete','type',cortype); %,'type','Spearman'
%[r p]=corr(Final{:,vars(1)},Final{:,vars(2)},'rows', 'pairwise','type',cortype); %,'type','Spearman'
r=cormat.R(vars(1),vars(2)); p=cormat.P(vars(1),vars(2));
%if isnan(r) | isnan(p)
[f, ~]=fit(corrs(Nn,vars(1)),corrs(Nn,vars(2)),'poly1');

NnA= isfinite(corrsADHD(:,vars(1))) & isfinite(corrsADHD(:,vars(2)));
%[rADHD pADHD]=corr(corrsADHD(NnA,vars(1)),corrsADHD(NnA,vars(2)),'rows', 'complete','type',cortype);
rADHD=cormatADHD.R(vars(1),vars(2)); pADHD=cormatADHD.P(vars(1),vars(2));
[fADHD, gofADHD]=fit(corrsADHD(NnA,vars(1)),corrsADHD(NnA,vars(2)),'poly1');
NnH= isfinite(corrsHealthy(:,vars(1))) & isfinite(corrsHealthy(:,vars(2)));
%[rHealthy pHealthy]=corr(corrsHealthy(NnH,vars(1)),corrsHealthy(NnH,vars(2)),'rows', 'complete','type',cortype);
rHealthy=cormatHealthy.R(vars(1),vars(2)); pHealthy=cormatHealthy.P(vars(1),vars(2));
[fHealthy, gofHealthy]=fit(corrsHealthy(NnH,vars(1)),corrsHealthy(NnH,vars(2)),'poly1');
ff=figure('Units', 'pixels', 'Position', [100, 100, 400, 400]); hold on;
xlim([min(corrs(Nn,vars(1)))-7.*ste(corrs(Nn,vars(1))) max(corrs(Nn,vars(1)))+7.*ste(corrs(Nn,vars(1)))])
ylim([min(corrs(Nn,vars(2)))-7.*ste(corrs(Nn,vars(2))) max(corrs(Nn,vars(2)))+7.*ste(corrs(Nn,vars(2)))])
% xline=[nanmin(corrs(Nn,vars(1))):0.1:nanmax(corrs(Nn,vars(1)))]';
% goff = predint(f,xline,0.95,'observation','off');
% gofADHD = predint(fADHD,xline,0.95,'observation','off');
% gofHealthy= predint(f,xline,0.95,'observation','off');
%ci = confint(f,0.95);
%ciup = polyval( ci(2,:),xline);
%cilow = polyval(ci(1,:),xline);
%corrlineADHD=lsline(ax) ; set(corrlineADHD,'color', ADHDcolor,'LineWidth',1); %plot(fADHD,'--')
%corr_a_line=lsline;
%corr_h_line=lsline;
corrlineg=plot(f) ; 
set(corrlineg,'color',[0 0 0 0.6],'LineWidth',1);
%gol=plot(xline,goff,'--','color', [0 0 0 0.6], 'LineWidth', 0.5);
%corrlineHealthy=lsline(h); set(corrlineHealthy, 'color', Healthycolor,'LineWidth',1); %plot(fHealthy,'--')
corr_a_line=plot(fADHD);
set(corr_a_line,'color', [ADHDcolor 0.95],'LineWidth',2,'LineStyle','--');
%golA=plot(xline,gofADHD,'--','color', [ADHDcolor 0.6], 'LineWidth', 0.5);
corr_h_line=plot(fHealthy);
set(corr_h_line, 'color', [Healthycolor 0.95],'LineWidth',2,'LineStyle','--');
%goH=plot(xline,gofHealthy,'--','color', [Healthycolor 0.6], 'LineWidth', 0.5);

a=scatter(corrsADHD(NnA,vars(1)),corrsADHD(NnA,vars(2)),120,'o','LineWidth',0.1,'MarkerFaceColor',ADHDcolor,'Markerfacealpha',.55,'MarkerEdgeColor',[ADHDcolor.*0.8]);
h=scatter(corrsHealthy(NnH,vars(1)),corrsHealthy(NnH,vars(2)),120,'o','LineWidth',0.1,'MarkerFaceColor',Healthycolor,'Markerfacealpha',.55,'MarkerEdgeColor',[Healthycolor.*0.8]);
title ({[cortype ' Correlation' newline strjoin({'General: R=' num2str(r,2), ' R^{2}' '=' num2str(r^2,2),'     p=' num2str(p,4)})...
    newline strjoin({[expgroup ': R= ' num2str(rADHD,2), '  R^{2}' '= ' num2str(rADHD^2,2),'      p= ' num2str(pADHD,4)]})...
    newline strjoin({[cntrgroup ': R= ' num2str(rHealthy,2),'  R^{2}' '= ' num2str(rHealthy^2,2),'      p= ' num2str(pHealthy,4)]})]});


xlabel(strrep(['Var ' num2str(corrvars(1)) '- ' Final.Properties.VariableNames{corrvars(1)}],'_',' '));
ylabel(strrep(['Var ' num2str(corrvars(2)) '- ' Final.Properties.VariableNames{corrvars(2)}],'_',' '));
set(gca,'FontName','Helvetica Neue','FontSize',10);
legend([a corr_a_line h corr_h_line corrlineg], [expgroup ' n=' num2str(size(corrsADHD(NnA,:),1))], [expgroup ' trend'], [cntrgroup ' n=' num2str(size(corrsHealthy(NnH,:),1))],[cntrgroup ' trend'], 'General trendline','FontSize',6)
set(ff,'Name',[cortype '_Corr_' expgroup '_' cntrgroup '_' Final.Properties.VariableNames{corrvars(1)} ' X ' Final.Properties.VariableNames{corrvars(2)}]);
warning('on')
%[r p]=corr(Final{:,vars(1)},Final{:,vars(2)},'rows', 'pairwise','type',cortype); %,'type','Spearman'
%[r p]=corr(Final{:,10},Final{:,30},'rows', 'complete') %,'type','Spearman'
%figure; scatter(Final{:,17},Final{:,30})



%cntrgroup, expgroup,
%    
%  ya=stats{stats.group==Groups(1),3};
%  ya=ya(~isnan(ya));
%  yb=stats{stats.group==Groups(2),3};
%  yb=yb(~isnan(yb));
% 
% colorsg=[ADHDcolor ; Healthycolor] ;
% 
% y=[ya;yb];
% 
% x=[cellstr(repmat([strrep(cntrgroup,'_',' ') '-' strrep(expgroup,'_',' ')],size(y,1),1))];
% group=[cellstr(repmat(strrep(Groups{1},'_',' '),size(ya,1),1)) ; cellstr(repmat(strrep(Groups{2},'_',' '),size(yb,1),1))];
% clear g
% g=gramm('x',x,'y',y,'color',group,'group',group);%,'color',bo,
% g.stat_violin('fill','transparent' ,'width',1 ,'dodge',0.5 );%'normalization','count', 'npoints',15,
% %g.stat_boxplot('width', 0.6 )
% g.stat_summary('type','sem' ,'geom','lines','dodge', 0.5 );%'quartile' '95percentile' 'std' 'fitnormalci' 
% %g.geom_jitter('width', 0.2,'dodge', 0.4 )
% %g.geom_vline('xintercept',1) 
% g.set_color_options('map',colorsg,'n_color',2,'n_lightness',1);
% g.set_text_options('font','Helvetica Neue','interpreter', 'tex',...
%     'legend_title_scaling', 1.1 , 'big_title_scaling', 1.2 );
% g.set_names('x','','y',['Amplitude (' '\mu' 'v)']);
% g.set_title(['p=' num2str(P) newline strrep(strjoin(stats.Properties.VariableNames(3)),'_',' ')]);
% figure('Position',[100 100 280 450]);
% %g.coord_flip();
% %g.set_order_options('color',1);
% 
% g.draw();

end
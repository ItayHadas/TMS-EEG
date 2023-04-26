function [ corMatF] = Corr_simp(Final, vars,alpha,cortype, disregard_zeros,disregard_num, rem_outlie) 
% corrvars=[3 73];[5 22];
% alpha=0.05
% cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
% Final=MST_FULL;
% Final.AgeAtTxStart=double(Final.AgeAtTxStart)
% Final.J_ROI63_64_time80_120 = cell2mat(Final.J_ROI63_64_time80_120)
%alpha=0.05

% vars=[56 70];
% alpha=0.05;
% cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
% Final=statsc; corrvars=vars;
% [correlations]=Corr_simp(Final, vars,alpha,cortype);



warning('off')
%vars=corrvars-2;
%vars=corrvars;
%groups=Final.group;

if sum(strcmpi(Final.Properties.VariableNames,'filt'))>0
    Final=Final(Final.filt==1,:);
end
tabdbl=~strcmp(varfun(@class,Final,'OutputFormat','cell'),'double');
fin2=table2cell(Final); fin2(:,tabdbl)=num2cell(nan(size(Final,1), sum(tabdbl)));
%corrs=Final{:,3:end} ;
corrs=cell2mat(fin2);
%corrs=table2array(Final(:,3:end));
%corrs(corrs==0)=nan;
%zeros removal
if strcmpi(disregard_zeros,'on')
%     corrsADHD(corrsADHD==0)=nan;
%     corrsHealthy(corrsHealthy==0)=nan;
    corrs(corrs==0)=nan;
end

if disregard_num~=0
    %corrsADHD(corrsADHD==disregard_num)=nan;
    %corrsHealthy(corrsHealthy==disregard_num)=nan;
    corrs(corrs==disregard_num)=nan;
end

%outlier removal
if strcmpi(rem_outlie,'on')
%     corrsADHD(isoutlier(corrsADHD))=nan;
%     corrsHealthy(isoutlier(corrsHealthy))=nan;
    corrs(isoutlier(corrs))=nan; %'mean' | 'quartiles' | 'grubbs' | 'gesd'
end


[cormat.R, cormat.P]=corr(corrs,'rows','pairwise', 'type', cortype);
%[cormat.R cormat.P]=corrcoef(corrs,'rows','pairwise');%'rows','pairwise'
%corrplot(X)
corMat=nan(size(cormat.R));
corMat(cormat.P<alpha)=cormat.R(cormat.P<alpha);
corMatF=array2table(corMat,'VariableNames',Final.Properties.VariableNames(1:end),'RowNames',Final.Properties.VariableNames(1:end)');
%corMatF=corMatF(sum(~isnan(corMatF{:,:}),2)>=1,sum(~isnan(corMatF{:,:}))>=1);
if size(vars,2)>1
%Healthycolor=[43/255 87/255 151/255];% 'b' ;%[166/255 166/255 166/255];
%co=brewermap(3,'RdBu');
%ADHDcolor=co(1,:);%'r'; %[0 0 0];
%Healthycolor=co(3,:);

co=brewermap(4,'PuOr'); %PuOr RdBu
%ADHDcolor=co(2,:);
%ADHDcolor1=co(1,:);%'r'; %[0 0 0];
Healthycolor=co(3,:); % 'b' ;%[166/255 166/255 166/255];
Healthycolor1=co(4,:);


%Nn= ~isnan(corrs(:,vars(1))) & ~isnan(corrs(:,vars(2)));
Nn= isfinite(corrs(:,vars(1))) & isfinite(corrs(:,vars(2)));

%[r, p]=corr(corrs(:,vars(1)),corrs(:,vars(2)),'rows', 'pairwise','type',cortype); %,'type','Spearman'
r=cormat.R(vars(1),vars(2)); p=cormat.P(vars(1),vars(2));
%r=cormat.R(vars(1),vars(2)); p=cormat.P(vars(1),vars(2));
[f, ~]=fit(corrs(Nn,vars(1)),corrs(Nn,vars(2)),'poly1');
figure; hold on;

%xlim([min(corrs(Nn,vars(1))) max(corrs(Nn,vars(1)))])
%ylim([min(corrs(Nn,vars(2))) max(corrs(Nn,vars(2)))])
xlim([min(corrs(Nn,vars(1)))-3.*ste(corrs(Nn,vars(1))) max(corrs(Nn,vars(1)))+3.*ste(corrs(Nn,vars(1)))])
ylim([min(corrs(Nn,vars(2)))-3.*ste(corrs(Nn,vars(2))) max(corrs(Nn,vars(2)))+3.*ste(corrs(Nn,vars(2)))])

corrline=plot(f); set(corrline,'color','k','LineWidth',1);
h=scatter(corrs(Nn,vars(1)),corrs(Nn,vars(2)),85,'o','LineWidth',0.75,'MarkerFaceColor',Healthycolor,'Markerfacealpha',.75,'MarkerEdgeColor',Healthycolor1);

title ([cortype ' Correlation  ' newline strjoin({'R=' num2str(r,2), ' R^{2}' '=' num2str(r^2,2),'  p=' num2str(p,4) '    N=' num2str(sum(Nn)) })]);
xlabel(strrep(['Var ' num2str(vars(1)) '- ' Final.Properties.VariableNames{vars(1)}],'_',' '));
ylabel(strrep(['Var ' num2str(vars(2)) '- ' Final.Properties.VariableNames{vars(2)}],'_',' '));
%str=sprintf(['R=%1.2f  R^{2}' '=%1.2f  p=%1.3f  N=%1.0f'],r, r^2, p, sum(Nn)) ;
%T=annotation('textbox',[0.42 0.05 0.4 0.27],'String',str,'FitBoxToText','on');
%T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
%set(T, 'fontsize', 11,'FontName','Helvetica Neue');
%set(T, 'fontsize', 11,'FontName','Helvetica Neue', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
%str=['f(x) = ' num2str(f.p1) ' *x + ' num2str(f.p2)]
%t=annotation('textbox',[0.42 0.05 0.4 0.27],'String',str,'FitBoxToText','on');
%t.FontName='Helvetica Neue'; t.FontSize = 11; t.FontWeight='bold';
hold off
end
warning('on')
end
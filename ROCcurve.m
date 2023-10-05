function [CM] = ROCcurve(Final4,vars,resp,clsfun)

%clearvars -except Final3 Final4
%Final4=Final3;
%Final4=statsc;
%Final4=outerjoin(statsc,comp1.MST.HAMD(:,[1 5:9]),'MergeKeys',true,'Keys', [1]);
%load('D:\Google_Drive\PhD Zangen\ADHD\avi''s paper\BIG_FINAL_DATA_table.mat'); Final4=BIG_Final(:,[1:23 26:27 30:34 46:40 42:47]);
% Final4=BIG_Final;
% Final4=statsc;
%vars=[10 14];
%vars=[5];
%resp=2
%vars=[27];
%resp=25;%F
%vars=[48];
%resp=2;%Final4(:,108);
%clsfun='logistic'; %logistic   linear
filtnan=all(isfinite(Final4{:,vars}),2);
predictors=Final4(filtnan,[vars]);
response=categorical(Final4{filtnan,resp});

if islogical(Final4{filtnan,resp})
   cati(1)= categorical(true); cati(2)=categorical(false);
else
   cati=unique(response);
end
%cat1='responders';
%cat2='nonresponders';
cat1=cati(1); cat2=cati(2);
%response(response=='true')=cat1;
%response(response=='false')=cat2;
%Final4.Var3=categorical(cellstr(Final4.Var3));
%cati=categorical({cat1; cat2});

%clearvars -except vars resp clsfun predictors response cat1 cat2 cati Final4
if strcmp(clsfun,'linear')
    clear classificationDiscriminant validationPredictions validationScores
classificationDiscriminant = fitcdiscr(...
     predictors, ...
     response, ...
     'DiscrimType', 'linear', ...   % diaglinear  linear
     'Gamma', 0, ...
     'FillCoeffs', 'off', ...
     'ClassNames', cati);
    
%[validationPredictions, validationScores] = resubPredict(ClassifierObject);
[validationPredictions, validationScores,~] = predict(classificationDiscriminant, predictors);

elseif strcmp(clsfun,'logistic')
    clear concatenatedPredictorsAndResponse validationScores GeneralizedLinearModel
    concatenatedPredictorsAndResponse=[predictors, table(response==cat1)]; %      table(response)
    
    GeneralizedLinearModel = fitglm(...
    concatenatedPredictorsAndResponse, ...
    'Distribution', 'binomial', ...
    'link', 'logit');

%  [~, validationScores] = predict(GeneralizedLinearModel, predictors);
[validationScores , ~] = predict(GeneralizedLinearModel, predictors);
validationPredictions(validationScores(:,1)>=0.5)=cat1;
validationPredictions(validationScores(:,1)<0.5)=cat2;
validationPredictions=validationPredictions';
  
end

% ROC
[x,y,T,AUC,OPTROCPT,~,~] = perfcurve(response,validationScores(:,1),cat1,'ProcessNaN','ignore',...
    'Alpha',0.05,'BootType','norm','NBoot',100,'UseNearest','off');
    %[~, AccuracyLDA, Thr] = perfcurve(response, validationScores(:,1), cat1,'yCrit','accu');
X=x(:,1); Y=y(:,1); 
fig=figure; hold on; grid on; fig.PaperUnits='centimeters'; fig.Units='centimeters';
fig.PaperPositionMode = 'auto';
patch('XData',[X; 1],'YData',[Y; 0],...
    'FaceVertexAlphaData',0.6,'FaceAlpha',0.5, 'LineStyle','none','FaceColor',[0.796 0.886 0.945],'EdgeColor','none');
patch('XData',[0 ;1; 1 ;0],'YData',[0; 1 ;0 ;0],'FaceVertexAlphaData',0.6,'FaceAlpha',0.9, 'LineStyle','none',...
    'FaceColor',[0.796 0.886 0.955],'EdgeColor','none');
plot(X,Y) %plot([0 1],[0 1])
xlabel('False positive rate (1-Specificity)') ; ylabel('True positive rate (Sensitivity)')
if size(AUC,2)>1
title([join(['ROC for Classification by ' clsfun ' Regression']) ...
    join(strrep(string(predictors.Properties.VariableNames),'_',' '))  join(['AUC ' num2str(AUC(1)) '   CI: ' num2str(AUC([2 3]))])])
else
    title([join(['ROC for Classification by ' clsfun ' Regression'])  ...
    join(strrep(string(predictors.Properties.VariableNames),'_',' '))  join(['AUC ' num2str(AUC(1))])])
end  

thresh=T(OPTROCPT(2)==y(:,1) & OPTROCPT(1)==x(:,1));
%choosing optimal threshold
validationPredictions2=validationPredictions; validationPredictions2(validationScores(:,1)>=thresh)=cati(1);
validationPredictions2(validationScores(:,1)<thresh)=cati(2);
if size(validationPredictions2,1)==size(filtnan,1)
validationPredictions2=validationPredictions2(filtnan');
end
[C2,ord]= confusionmat(cellstr(response),cellstr(validationPredictions2),'Order',cellstr(cati)); 
TP=C2(1,1); TN=C2(2,2); FP=C2(2,1); FN=C2(1,2); accuracy=(TP+TN)./(TP+TN+FP+FN);
sensitivity=TP./(TP+FN); specificity=TN./(TN+FP); Likelihood_pos=sensitivity./(1-specificity);
Likelihood_neg=(1-sensitivity)./specificity;

plot([1-specificity],[sensitivity],'-s','MarkerSize',9,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
str={['Accuracy = ' num2str(accuracy)], ['Sensitivity = ' num2str(sensitivity)],...
    ['Specificity = ' num2str(specificity)],['Likelihood ratio positive = ' num2str(Likelihood_pos)],...
    ['Likelihood ratio negative = ' num2str(Likelihood_neg)]}; 
t=annotation('textbox',[0.45 0.08 0.4 0.29],'String',str,'FitBoxToText','on');
t.FontName='Helvetica Neue'; t.FontSize = 11; t.FontWeight='bold'; t.LineStyle='none';
hold off
figure('position', [500 200 520 210])
CM=confusionchart(cellstr(response),cellstr(validationPredictions2));
sortClasses(CM,ord)
if size(AUC,2)>1
title([join(['Classification by ' clsfun ' Regression']) ...
    join(strrep(string(predictors.Properties.VariableNames),'_',' '))  join(['AUC ' num2str(AUC(1)) '   CI: ' num2str(AUC([2 3]))])])
else
    title([join(['Classification by ' clsfun ' Regression'])  ...
    join(strrep(string(predictors.Properties.VariableNames),'_',' '))  join(['AUC ' num2str(AUC(1))])])
end  
end


%% junk

% confussion matrix for the optimal ROC point
%[C1,~]= confusionmat(cellstr(response),cellstr(validationPredictions),'Order',cellstr(cati)); %,'Order',string(cati)
%cp = classperf(cellstr(response));
%TP1=C1(1,1); TN1=C1(2,2); FP1=C1(2,1); FN1=C1(1,2); %accuracy1 = (TP1+TN1)./(TP1+TN1+FP1+FN1); sensitivity1=TP1./(TP1+FN1); specificity1=TN1./(TN1+FP1);
%plot([1-specificity1],[sensitivity1],'-s','MarkerSize',9,'MarkerEdgeColor','yellow','MarkerFaceColor',[.5 .9 .9])
%plot(OPTROCPT(1),OPTROCPT(2),'r*')
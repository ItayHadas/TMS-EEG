

%% rearanging 

a=readtable('D:\WORKINGset-D\STAR_d_suicidality\1_DEMOGRAPHICSdm011.xlsx');
b=readtable('D:\WORKINGset-D\STAR_d_suicidality\2_qids012.xlsx');
b(strcmpi(b.FormUsed_assessmentName,''),:)=[]; b(strcmpi(b.TreatmentLevel, 'Level 2.1'),:)=[];  b.Properties.VariableNames{1}='SubjectID';
c=readtable('D:\WORKINGset-D\STAR_d_suicidality\3_MEDHX_ccv011.xlsx');
d=readtable('D:\WORKINGset-D\STAR_d_suicidality\4_hrsd011.xlsx');




unique(b.TreatmentLevel)
size(unique(a{:,1}))
b(strcmpi(b{:,5},'Level 3'),[4 6 7])
a=sortrows(a,1);

c(strcmpi(c.TreatmentLevel,'Follow-Up'),:)
unstack(a,'TreatmentLevel','AgeInMonthsAtTheTimeOfTheInterview_test_sampling_imaging_','AggregationFunction',@isnan)

{'Enrollment'}
{'Follow up' }
{'Level 1'   }
{'Level 2'   }
{'Level 2A'  }
{'Level 3'   }
{'Level 4'   }


%% QIDS re-arangement 
%b=readtable('D:\WORKINGset-D\STAR_d_suicidality\2_qids012.xlsx');
%b=readtable('D:\WORKINGset-D\STAR_d_suicidality\qids01.xlsx');
b.Properties.VariableNames{5} = 'TreatmentLevel'; b.Properties.VariableNames{1}='SubjectID';
b(strcmpi(b.TreatmentLevel, 'Level 1'),:).TreatmentLevel=repmat({'Level1'},size(b(strcmpi(b.TreatmentLevel, 'Level 1'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 2.1'),:).TreatmentLevel=repmat({'Level2A'},size(b(strcmpi(b.TreatmentLevel, 'Level 2.1'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level_2-1'),:).TreatmentLevel=repmat({'Level2A'},size(b(strcmpi(b.TreatmentLevel, 'Level_2-1'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 2'),:).TreatmentLevel=repmat({'Level2'},size(b(strcmpi(b.TreatmentLevel, 'Level 2'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 3'),:).TreatmentLevel=repmat({'Level3'},size(b(strcmpi(b.TreatmentLevel, 'Level 3'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 4'),:).TreatmentLevel=repmat({'Level4'},size(b(strcmpi(b.TreatmentLevel, 'Level 4'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Follow-Up'),:).TreatmentLevel=repmat({'zFollow-Up'},size(b(strcmpi(b.TreatmentLevel, 'Follow-Up'),:),1),1);
b(strcmp(b.FormUsed_assessmentName,''),:)=[]; b(strcmp(b.TreatmentLevel,''),:)=[];
b=sortrows(b,4); b=sortrows(b,5); b=sortrows(b,1);

% histogram(b.DaysSinceBaseline); hold on;
% export_fig sample_weeks.pdf -q101 -painters -append
%jjj=b ;jjj(strcmpi(jjj.TreatmentLevel, 'zFollow-Up'),:)=[];
% histogram(floor(jjj.DaysSinceBaseline/7))
%b_backup=b;

% Take the highest score of items 1-4, 5
% the highest score of items 6-9, 10 11 12 13 14
% the highest score of items 15-16,
% and then all the remaining scores
% (so overall you are adding up 10 numbers to get the total score for QIDS).
% re-score the total QIDS score
for i=1:size(b,1)%[unique(b.SubjectID)]' %i=2846
    %if isnan(b.QIDSTotalScore(i))
    b.QIDSTotalScore2(i)=max(b{i,6:9})+max(b{i,11:14})+max(b{i,20:21})+b{i,10}+nansum(b{i,15:19});
    %end
end

%% histogram for a cetain condition in the table
for i=[unique(b.SubjectID)]'
  iii(i)=std(table2array(demog([demog.src_subject_id==i],[6])));
end
hist(iii)

%% meds re-arangement 
% unique(m.TreatmentLevel)
m(strcmpi(m.TreatmentLevel, 'Level 1'),:).TreatmentLevel=repmat({'Level1'},size(m(strcmpi(m.TreatmentLevel, 'Level 1'),:),1),1);
m(strcmpi(m.TreatmentLevel, 'Level 2'),:).TreatmentLevel=repmat({'Level2'},size(m(strcmpi(m.TreatmentLevel, 'Level 2'),:),1),1);
m(strcmpi(m.TreatmentLevel, 'Level 2A'),:).TreatmentLevel=repmat({'Level2A'},size(m(strcmpi(m.TreatmentLevel, 'Level 2A'),:),1),1);
m(strcmpi(m.TreatmentLevel, 'Level 3'),:).TreatmentLevel=repmat({'Level3'},size(m(strcmpi(m.TreatmentLevel, 'Level 3'),:),1),1);
m(strcmpi(m.TreatmentLevel, 'Level 4'),:).TreatmentLevel=repmat({'Level4'},size(m(strcmpi(m.TreatmentLevel, 'Level 4'),:),1),1);
m.week(m.week==0.1,:)=0;
%unique(m.week)
warning off;

%% Cognitive Meds re-arangement 

cogmed(strcmp(cogmed.level,''),:)=[];
cogmed(strcmpi(cogmed.level, 'Level 1'),:).level=repmat({'Level1'},size(cogmed(strcmpi(cogmed.level, 'Level 1'),:),1),1);
cogmed(strcmpi(cogmed.level, 'Level 2.1'),:).level=repmat({'Level2A'},size(cogmed(strcmpi(cogmed.level, 'Level 2.1'),:),1),1);
cogmed(strcmpi(cogmed.level, 'Level_2-1'),:).level=repmat({'Level2A'},size(cogmed(strcmpi(cogmed.level, 'Level_2-1'),:),1),1);
cogmed(strcmpi(cogmed.level, 'Level 2'),:).level=repmat({'Level2'},size(cogmed(strcmpi(cogmed.level, 'Level 2'),:),1),1);
cogmed(strcmpi(cogmed.level, 'Level 3'),:).level=repmat({'Level3'},size(cogmed(strcmpi(cogmed.level, 'Level 3'),:),1),1);
cogmed(strcmpi(cogmed.level, 'Level 4'),:).level=repmat({'Level4'},size(cogmed(strcmpi(cogmed.level, 'Level 4'),:),1),1);
cogmed(strcmpi(cogmed.level, 'Follow-Up'),:).level=repmat({'zFollow-Up'},size(cogmed(strcmpi(cogmed.level, 'Follow-Up'),:),1),1);
cogmed=sortrows(cogmed,4); cogmed=sortrows(cogmed,7); 


writetable(cogmed,'D:\WORKINGset-D\STAR_d_suicidality\cognitive_tx2.xls')

%% COMPLETE data integration

%writetable(b,'D:\WORKINGset-D\STAR_d_suicidality\qids02.xlsx');
%writetable(m,'D:\WORKINGset-D\STAR_d_suicidality\STARD-medications.xlsx');

b=readtable('D:\WORKINGset-D\STAR_d_suicidality\qids02.xlsx');
d=readtable('D:\WORKINGset-D\STAR_d_suicidality\hrsd02.xlsx');
m=readtable('D:\WORKINGset-D\STAR_d_suicidality\STARD-medications.csv','Format','%f %s %f %s %s %s %s %f %f %f %f');
demog=readtable('D:\WORKINGset-D\STAR_d_suicidality\DEMOGRAPHICS_for_itay.xlsx');
hst=readtable('D:\WORKINGset-D\STAR_d_suicidality\clinical_history.xls');
ppsy=readtable('D:\WORKINGset-D\STAR_d_suicidality\past_psychotropic.xls');
raceinf=readtable('D:\WORKINGset-D\STAR_d_suicidality\Race_info.xls');
suici=readtable('D:\WORKINGset-D\STAR_d_suicidality\suicide_data.xls');
suici.Properties.VariableNames{27} = 'Attempted_suicide';
cogmed=readtable('D:\WORKINGset-D\STAR_d_suicidality\cognitive_tx2.xls');


temp=table; temp.SubjectID=unique(b.SubjectID); %temp=[temp array2table(nan(size(temp,1),11))];
clear i
for i=[unique(b.SubjectID)]'
    
    %% demographics
    pivD=table;
    pivD=demog([demog.src_subject_id==i],[4 6 7 15 20 21 23 26 27 28 41]);
    if size(pivD,1)>1
        ll=[((min(sum(isnan(pivD{:,[4 5 7 8 9 10]})')'))==(sum(isnan(pivD{:,[4 5 7 8 9 10]})')'))==0];
        income=nanmean(pivD{:,11});
        age=nanmean(pivD{:,2});
        pivD(ll,:)=[];
        if size(pivD,1)>1; pivD=pivD(end,:); end
        pivD(:,11)={round(income)}; pivD(:,2)={round(age)};
    end
    if isempty(pivD)
        tt=cell2table({nan,nan,'',nan,nan,'',nan,nan,nan,nan,nan});%({[],[],'',[],[],'',[],[],[],[],[]});
        tt.Properties.VariableNames={'src_subject_id','interview_age','gender','marital','educat','grade_highed','empl','publica','medicaid','privins','totincom'};
        pivD=tt;
        clear tt
    end
    if i==1
        pivD=[pivD ; repmat({nan,nan,'',nan,nan,'',nan,nan,nan,nan,nan},size(temp,1)-1,1)];%[pivD ; cell(size(temp,1)-1,11)];
        temp=[temp pivD(:,2:end)];
        % clear pivD
    else
        temp([temp.SubjectID==i],[2:11])=pivD(:,2:end);
    end
    
    %% clinical history
    pivhs=table;
    pivhs=hst([hst.src_subject_id==i],[4 11:38 42 46 50 54 61 110]);
    if size(pivhs,1)>1
        ll=[((min(sum(isnan(pivhs{:,[2:35]})')'))==(sum(isnan(pivhs{:,[2:35]})')'))==0];
        pivhs(ll,:)=[];
        if size(pivhs,1)>1; pivhs=pivhs(end,:); end
    end
    if isempty(pivhs)
        pivhs=array2table(nan(1,size(pivhs,2)),'VariableNames',{'src_subject_id','dage','epino','episode_date','ai_none','ai_def','ai_na','alcoh','amphet','cannibis','opioid','pd_ag','pd_noag','specphob','soc_phob','ocd_phx','psd','gad_phx','axi_oth','aii_none','aii_def','aii_na','pd_border','pd_depend','pd_antis','pd_paran','pd_nos','axii_oth','dep','bip','alcohol','drug_phx','suic_phx','ax_cocaine','etcmats'});
    end
    if i==1
        pivhs=[pivhs ; array2table(nan(size(temp,1)-1,size(pivhs,2)),'VariableNames',pivhs.Properties.VariableNames)];
        temp=[temp pivhs(:,2:end)];
    else
        temp([temp.SubjectID==i],[12:45])=pivhs(:,2:end);
    end
    
    %% past_psychotropic
    pivpsy=table;
    pivpsy=ppsy([ppsy.src_subject_id==i],[4 10]);
    if size(pivpsy,1)>1
        ll=[((min(sum(isnan(pivpsy{:,[2]})')'))==(sum(isnan(pivpsy{:,[2]})')'))==0];
        pivpsy(ll,:)=[];
        if size(pivpsy,1)>1; pivpsy=pivpsy(end,:); end
    end
    if isempty(pivpsy)
        pivpsy=array2table(nan(1,size(pivpsy,2)),'VariableNames',{'src_subject_id','psmed'});
    end
    if i==1
        pivpsy=[pivpsy ; array2table(nan(size(temp,1)-1,size(pivpsy,2)),'VariableNames',pivpsy.Properties.VariableNames)];
        temp=[temp pivpsy(:,2:end)];
    else
        temp([temp.SubjectID==i],[46])=pivpsy(:,2:end);
    end
    
    %% race info
    pivrace=table;
    pivrace=raceinf([raceinf.src_subject_id==i],[4 11]);
    if size(pivrace,1)>1
        ll=[((min(sum(isnan(pivrace{:,[2]})')'))==(sum(isnan(pivrace{:,[2]})')'))==0];
        pivrace(ll,:)=[];
        if size(pivrace,1)>1; pivrace=pivrace(end,:); end
    end
    if isempty(pivrace)
        pivrace=array2table(cell(1,size(pivrace,2)),'VariableNames',{'src_subject_id','race'});
    end
    if i==1
        pivrace=[pivrace ; array2table(cell(size(temp,1)-1,size(pivrace,2)),'VariableNames',pivrace.Properties.VariableNames)];
        temp=[temp pivrace(:,2:end)];
    else
        temp([temp.SubjectID==i],[47])=pivrace(:,2:end);
    end
    
    %% Suicide data suici (only level 1)
    pivsui=table;
    pivsui=sortrows(suici(suici.src_subject_id==i,[5 27]),1);
    
    if size(pivsui,1)>1
        ll=[((min(sum(isnan(pivsui{:,[2]})')'))==(sum(isnan(pivsui{:,[2]})')'))==0];
        pivsui(ll,:)=[];
        if size(pivsui,1)>1; pivsui=pivsui(end,:); end
    end
    if isempty(pivsui)
        pivsui=array2table(nan(1,size(pivsui,2)),'VariableNames',{'src_subject_id','Attempted_suicide'});
    end
    if i==1
        pivsui=[pivsui ; array2table(nan(size(temp,1)-1,size(pivsui,2)),'VariableNames',pivsui.Properties.VariableNames)];
        temp=[temp pivsui(:,2:end)];
    else
        temp([temp.SubjectID==i],[48])=pivsui(:,2);
    end
    
    %%
    Treat= unique(b(b.SubjectID==i,:).TreatmentLevel);
    for k=[1:size(Treat,1)]
        Day= sort(unique(b([b.SubjectID==i & strcmpi(b.TreatmentLevel,Treat(k))],:).DaysSinceBaseline));
        days1=[1:size(Day,1)];
        
        for j=days1   %[1:size(Day,1)]
    %% QIDS
            form=unique(b([b.SubjectID==i & strcmpi(b.TreatmentLevel,Treat(k)) & b.DaysSinceBaseline==Day(j) ],:).FormUsed_assessmentName);
            for l=[1:size(form,1)] %Clinician/Self Forms
            
                piv=table;
                piv=b([b.SubjectID==i & strcmpi(b.TreatmentLevel,Treat(k)) & b.DaysSinceBaseline==Day(j) & strcmpi(b.FormUsed_assessmentName,form(l))],: );
                
                if size(piv,1)>1
                    ll=[((min(sum(isnan(piv{:,[17 33]})')'))==(sum(isnan(piv{:,[17 33]})')'))==0];
                    piv(ll,:)=[];
                    if size(piv,1)>1; piv=piv(end,:); end
                end
                piv=piv(:,[1 3 4 5 17 33]);
                Week=sort(floor(Day/7)-floor(Day(1)/7));
                if Week(1)==0
                    weeknum=[0:size(Week,1)-1];
                else
                    weeknum=[1:size(Week,1)]; % weeknum=1:size(Week,1);
                end
                
                var1= strrep(strrep(['QIDSSuIdea_' char(Treat(k)) '_week' num2str(weeknum(j)) '_' form{[l]}(1:4)],' ','_'),'-','_');
                var2= strrep(strrep(['QIDS_sum_' char(Treat(k)) '_week' num2str(weeknum(j)) '_' form{[l]}(1:4)],' ','_'),'-','_');
                if sum(strcmpi(temp.Properties.VariableNames,var1))==0
                    temp.temp1=nan(size(temp,1),1); %table('size',[size(temp,1) 1],'VariableTypes',{'double'},'VariableNames',{'temp1'})
                    temp([temp.SubjectID==i],:).temp1=piv.QIDSSuicidalIdeation;
                    temp.Properties.VariableNames{'temp1'} = var1;
                else
                    temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,var1))={piv.QIDSSuicidalIdeation};
                end
                if sum(strcmpi(temp.Properties.VariableNames,var2))==0
                    temp.temp2=nan(size(temp,1),1);
                    temp([temp.SubjectID==i],:).temp2=piv.QIDSTotalScore2;
                    temp.Properties.VariableNames{'temp2'} = var2;
                else
                    temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,var2))={piv.QIDSTotalScore2};
                end
                
            end %Clinician/Self Forms
          
    %% Medications
                      
            med_wk=sort(m{[m.SubjectID==i & strcmpi(m.TreatmentLevel,Treat(k)) ],3},'ascend');
            if ~isempty(med_wk)
                mesure_dist=pdist2(med_wk,Week);
                if j==days1(end)
                    gap1=4.1;
                elseif j==days1(1) && Week(j)==0 && med_wk(1)==0
                    gap1=0.1;
                else
                    gap1=1.1;
                end
                
                if min(mesure_dist(:,j))<gap1
                    [yo, II]=min(mesure_dist(:,j));
                    pivm=table;
                    pivm=sortrows([ m([m.SubjectID==i & strcmpi(m.TreatmentLevel,Treat(k)) ],[3:end])],1);
                    pivm=pivm(II,:);
                
            
            
            if size(pivm,1)>1
                [~,loc5]=min(sum(cellfun(@isempty,(table2cell(pivm(:,[1 5])))),2));
                pivm=pivm(loc5,:);
                if size(pivm,1)>1; pivm=pivm(1,:); end
            end
            if ~isempty(pivm)
                pivm=[piv(:,3) pivm];
                
                var3=['_' char(Treat(k)) '_week' num2str(weeknum(j)) '_'];
                pivm.Properties.VariableNames=strrep(strcat(pivm.Properties.VariableNames,var3),'-','_');
                if sum(contains(temp.Properties.VariableNames,pivm.Properties.VariableNames(1)))==0
                    temp3=num2cell(nan(size(temp,1),size(pivm,2))); %array2table(nan(size(temp,1),size(pivm,2)));
                    temp3([temp.SubjectID==i],:)=table2cell(pivm);
                    temp3=array2table(temp3);
                    temp3.Properties.VariableNames=pivm.Properties.VariableNames;
                    temp=[temp temp3];
                    
                else
                    temp{[temp.SubjectID==i],pivm.Properties.VariableNames}=table2cell(pivm);
                end
            end    
            end
            end
    %% Cognitive Medicine
            pivcog=table;
            pivcog=sortrows([ cogmed([cogmed.src_subject_id==i & strcmpi(cogmed.level,Treat(k)) ],[4 7 8 11 12])],1);
            if ~isempty(pivcog)
            if size(pivcog,1)>1
                ll=[((min(sum(isnan(pivcog{:,[4 5]})')'))==(sum(isnan(pivcog{:,[4 5]})')'))==0];
                pivcog(ll,:)=[];
                if size(pivcog,1)>1
                    mesure_dist=pdist2(pivcog.days_baseline,Day);
                    [~,i1]= min(mesure_dist,[],1);
                    pivcog=pivcog(i1(1),:);
                end
            end
            [~,i2]= min(mesure_dist,[],2);
            if j==i2(1)
                pivcog.Properties.VariableNames{4} = strrep(strrep(['CogSwitch_' char(Treat(k)) '_week' num2str(weeknum(j))],' ','_'),'-','_');
                pivcog.Properties.VariableNames{5} = strrep(strrep(['CogAugmentation_' char(Treat(k)) '_week' num2str(weeknum(j))],' ','_'),'-','_');
                
                if sum(strcmpi(temp.Properties.VariableNames,pivcog.Properties.VariableNames{4}))==0
                    temp.temp1=nan(size(temp,1),1);  temp.temp2=nan(size(temp,1),1);
                    temp([temp.SubjectID==i],:).temp1=pivcog{:,[4]};
                    temp([temp.SubjectID==i],:).temp2=pivcog{:,[5]};
                    temp.Properties.VariableNames{'temp1'} = pivcog.Properties.VariableNames{4};
                    temp.Properties.VariableNames{'temp2'} = pivcog.Properties.VariableNames{5};
                else
                    temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,pivcog.Properties.VariableNames{4}))={pivcog{:,[4]}};
                    temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,pivcog.Properties.VariableNames{5}))={pivcog{:,[5]}};
                end
            end
            end
        end
    end
    %% HRSD
    Treat= unique(d(d.SubjectID==i,:).TreatmentLevel);
    for k=[1:size(Treat,1)]
         DayH= sort(unique(d([d.SubjectID==i & strcmpi(d.TreatmentLevel,Treat(k))],:).DaysSinceBaseline));
         daysH=[1:size(DayH,1)];
         if ~isempty(daysH)
         for j=daysH
             piv=table;
             piv=d([d.SubjectID==i & strcmpi(d.TreatmentLevel,Treat(k)) & d.DaysSinceBaseline==DayH(j)],[1 3 4:22] );
             if size(piv,1)>1
                 ll=[((min(sum(isnan(piv{:,[3 4]})')'))==(sum(isnan(piv{:,[3 4]})')'))==0];
                 if sum(ll)<1
                     ll(1:end-1,:)=1;
                 end
                 piv(ll,:)=[];
                 if size(piv,1)>1; piv=piv(end,:); end
             end
             Week=sort(floor(DayH/7)-floor(DayH(1)/7));
                if Week(1)==0
                    weeknum=[0:size(Week,1)-1];
                else
                    weeknum=[1:size(Week,1)]; % weeknum=1:size(Week,1);
                end
            var1= strrep(strrep(['HRSSuicide_' char(Treat(k)) '_week' num2str(weeknum(j)) '_'  ],' ','_'),'-','_');
            var2= strrep(strrep(['HRSTotalScore_' char(Treat(k)) '_week' num2str(weeknum(j)) '_' ],' ','_'),'-','_');
            
            if sum(strcmpi(temp.Properties.VariableNames,var1))==0
                temp.temp1=nan(size(temp,1),1);
                temp([temp.SubjectID==i],:).temp1=cell2mat({piv.HRSSuicide});
                temp.Properties.VariableNames{'temp1'} = var1;
            else
                temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,var1))={piv.HRSSuicide};
            end
            if sum(strcmpi(temp.Properties.VariableNames,var2))==0
                temp.temp2=nan(size(temp,1),1);
                temp([temp.SubjectID==i],:).temp2=cell2mat({piv.HRSTotalScore_recorded_});
                temp.Properties.VariableNames{'temp2'} = var2;
            else
                temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,var2))={piv.HRSTotalScore_recorded_};
            end
            
         end
         end
    end
end


writetable(temp,['D:\WORKINGset-D\STAR_d_suicidality\QIDStemp_12.xlsx']);
%% SORT
% temp=readtable('D:\WORKINGset-D\STAR_d_suicidality\QIDStemp6.xlsx');
%warning on;
tempbak=temp;
% temp=tempbak;
temp=sortrows(temp,1);
temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2_'))=strrep(temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2_')),'Level2_','Level2a_');
temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2A_'))=strrep(temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2A_')),'Level2A_','Level2B_');
temp=[temp(:,[1]) temp(:,sort_nat(temp.Properties.VariableNames(2:end)))];

wk0=contains(temp.Properties.VariableNames(1:end),'week0_');
wk1=contains(temp.Properties.VariableNames(1:end),'week1_');
wk2=contains(temp.Properties.VariableNames(1:end),'week2_');
wk3=contains(temp.Properties.VariableNames(1:end),'week3_');
wk4=contains(temp.Properties.VariableNames(1:end),'week4_');
wk5=contains(temp.Properties.VariableNames(1:end),'week5_');
wk6=contains(temp.Properties.VariableNames(1:end),'week6_');
wk7=contains(temp.Properties.VariableNames(1:end),'week7_');
wk8=contains(temp.Properties.VariableNames(1:end),'week8_');
wk9=contains(temp.Properties.VariableNames(1:end),'week9_');
wk10=contains(temp.Properties.VariableNames(1:end),'week10_');
wk11=contains(temp.Properties.VariableNames(1:end),'week11_');
wk12=contains(temp.Properties.VariableNames(1:end),'week12_');
wk13=contains(temp.Properties.VariableNames(1:end),'week13_');
wk14=contains(temp.Properties.VariableNames(1:end),'week14_');
%sum(contains(temp.Properties.VariableNames(1:end),'week15_'));
temp=[temp(:,[1]) temp(:,[wk0]) temp(:,[wk1]) temp(:,[wk2]) temp(:,[wk3]) temp(:,[wk4]) temp(:,[wk5]) temp(:,[wk6]) temp(:,[wk7]) temp(:,[wk8]) temp(:,[wk9]) temp(:,[wk10]) temp(:,[wk11]) temp(:,[wk12]) temp(:,[wk13]) temp(:,[wk14])];



lvl1=contains(temp.Properties.VariableNames(1:end),'Level1_');
lvl2=contains(temp.Properties.VariableNames(1:end),'Level2_');
lvl2A=contains(temp.Properties.VariableNames(1:end),'Level2a_');
lvl2B=contains(temp.Properties.VariableNames(1:end),'Level2B_');
lvl3=contains(temp.Properties.VariableNames(1:end),'Level3_');
lvl4=contains(temp.Properties.VariableNames(1:end),'Level4_');
flup=contains(temp.Properties.VariableNames(1:end),'zFollow_Up_');
temp=[temp(:,[1]) temp(:,[lvl1]) temp(:,[lvl2]) temp(:,[lvl2A]) temp(:,[lvl2B]) temp(:,[lvl3]) temp(:,[lvl4]) temp(:,[flup])];





temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2a_'))=strrep(temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2a_')),'Level2a_','Level2_');
temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2B_'))=strrep(temp.Properties.VariableNames(contains(temp.Properties.VariableNames(1:end),'Level2B_')),'Level2B_','Level2A_');
QIDS=temp;
%med=temp(:,logical([ones(1,1) contains(temp.Properties.VariableNames(2:end),'medication')]));
%clin=temp(:,logical([ones(1,1) contains(temp.Properties.VariableNames(2:end),'_Clin')]));
%self=temp(:,logical([1 contains(temp.Properties.VariableNames(2:end),'_Self')]));
%QIDS=outerjoin(clin,self,'MergeKeys',true);
%QIDS=outerjoin(med,QIDS,'MergeKeys',true);
writetable(QIDS,'D:\WORKINGset-D\STAR_d_suicidality\QIDS2.xlsx')
clearvars Day form i j k l ll piv QIDS temp Treat var1 var2

%% General LOCF
% clear all

opts = detectImportOptions(['D:\WORKINGset-D\STAR_d_suicidality\QIDStemp_12.xlsx']);
varchange=opts.VariableNames([contains(opts.VariableNames,{'medication'}) & contains(opts.VariableNames,{'name'})]);
varchange2=opts.VariableNames([contains(opts.VariableNames,{'QIDS_sum'}) | contains(opts.VariableNames,{'QIDSSuIdea'}) | contains(opts.VariableNames,{'dosage'})]);
varchange3=opts.VariableNames([contains(opts.VariableNames,{'medication'}) & contains(opts.VariableNames,{'name'})]);
varchange4=opts.VariableNames([contains(opts.VariableNames,{'HRSSuicide'}) | contains(opts.VariableNames,{'HRSTotalScore'}) | contains(opts.VariableNames,{'dosage'})]);
varchange5=opts.VariableNames([contains(opts.VariableNames,{'CogSwitch'}) | contains(opts.VariableNames,{'CogAugmentation'})]);

opts = setvartype(opts,varchange,'string');
opts = setvartype(opts,varchange2,'double');
opts = setvartype(opts,varchange3,'string');
opts = setvartype(opts,varchange4,'double');
opts = setvartype(opts,varchange5,'double');

temp=readtable(['D:\WORKINGset-D\STAR_d_suicidality\QIDStemp_12.xlsx'],opts);
temp_bak=temp;

level={}; level(1,1)= {'Enrollment' }; level(2,1)= {'Level1_' }; level(3,1)={'Level2_'} ; level(4,1)={'Level2A_'} ; level(5,1)={'Level3_'} ; level(6,1)={'Level4_'} ; level(7,1)={'zFollow_Up_'};



levelind(1,:)=contains(temp.Properties.VariableNames,level(1,1));
levelind(2,:)=contains(temp.Properties.VariableNames,level(2,1));
levelind(3,:)=contains(temp.Properties.VariableNames,level(3,1));
levelind(4,:)=contains(temp.Properties.VariableNames,level(4,1));
levelind(5,:)=contains(temp.Properties.VariableNames,level(5,1));
levelind(6,:)=contains(temp.Properties.VariableNames,level(6,1));
levelind(7,:)=contains(temp.Properties.VariableNames,level(7,1));

qidssum=contains(temp.Properties.VariableNames,{'QIDS_sum'});
QIDSSuIdea=contains(temp.Properties.VariableNames,{'QIDSSuIdea'});
clin=contains(temp.Properties.VariableNames,{'_Clin'});
self=contains(temp.Properties.VariableNames,{'_Self'});
HRSsui=contains(temp.Properties.VariableNames,{'HRSSuicide'});
HRStot=contains(temp.Properties.VariableNames,{'HRSTotalScore'});
COGswch=contains(temp.Properties.VariableNames,{'CogSwitch'});
COGaug=contains(temp.Properties.VariableNames,{'CogAugmentation'});

med_dose1= [contains(temp.Properties.VariableNames,{'medication1'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose2= [contains(temp.Properties.VariableNames,{'medication2'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose3= [contains(temp.Properties.VariableNames,{'medication3'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose4= [contains(temp.Properties.VariableNames,{'medication4'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_name1= [contains(temp.Properties.VariableNames,{'medication1'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name2= [contains(temp.Properties.VariableNames,{'medication2'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name3= [contains(temp.Properties.VariableNames,{'medication3'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name4= [contains(temp.Properties.VariableNames,{'medication4'}) & contains(temp.Properties.VariableNames,{'name'})];



for i=1:size(temp,1)
    
    for k=1:size(levelind,1)
        
        clear locq
        locq{1,:}=temp{i,[levelind(k,:) & qidssum & clin]};
        %temp(i,[levelind(k,:) & qidssum & clin]);
        locq{2,:}=temp{i,[levelind(k,:) & qidssum & self]};
        %temp(i,[levelind(k,:) & qidssum & self]);
        locq{3,:}=temp{i,[levelind(k,:) & QIDSSuIdea & clin]};
        %temp(i,[levelind(k,:) & QIDSSuIdea & clin]);
        locq{4,:}=temp{i,[levelind(k,:) & QIDSSuIdea & self]};
        %temp(i,[levelind(k,:) & QIDSSuIdea & self]);
        locq{5,:}=temp{i,[levelind(k,:) & HRSsui]};
        locq{6,:}=temp{i,[levelind(k,:) & HRStot]};
        locq{7,:}=temp{i,[levelind(k,:) & COGswch]};
        locq{8,:}=temp{i,[levelind(k,:) & COGaug]};
        
        v=LOCF(locq);
        clear locq
        if ~isempty(v{1}), temp{i,[levelind(k,:) & qidssum & clin]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & qidssum & self]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & QIDSSuIdea & clin]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & QIDSSuIdea & self]}=v{4}; end
        if ~isempty(v{5}), temp{i,[levelind(k,:) & HRSsui]}=v{5}; end
        if ~isempty(v{6}), temp{i,[levelind(k,:) & HRStot]}=v{6}; end
        if ~isempty(v{7}), temp{i,[levelind(k,:) & COGswch]}=v{7}; end
        if ~isempty(v{8}), temp{i,[levelind(k,:) & COGaug]}=v{8}; end
        clear v
        %
        clear locd
        locd{1,:}=(temp{i,[levelind(k,:) & med_dose1]});
        locd{2,:}=(temp{i,[levelind(k,:) & med_dose2]});
        locd{3,:}=(temp{i,[levelind(k,:) & med_dose3]});
        locd{4,:}=(temp{i,[levelind(k,:) & med_dose4]});
        v=LOCF(locd);
        if ~isempty(v{1}), temp{i,[levelind(k,:) & med_dose1]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & med_dose2]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & med_dose3]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & med_dose4]}=v{4}; end
        clear v
        %
        clear locn
        locn{1,:}=temp{i,[levelind(k,:) & med_name1]};%table2cell(temp(i,[levelind(k,:) & med_name1]));
        locn{2,:}=temp{i,[levelind(k,:) & med_name2]};%table2cell(temp(i,[levelind(k,:) & med_name2]));
        locn{3,:}=temp{i,[levelind(k,:) & med_name3]};%table2cell(temp(i,[levelind(k,:) & med_name3]));
        locn{4,:}=temp{i,[levelind(k,:) & med_name4]};%table2cell(temp(i,[levelind(k,:) & med_name4]));
        v=LOCF(locn);
                
        if ~isempty(v{1}), temp{i,[levelind(k,:) & med_name1]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & med_name2]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & med_name3]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & med_name4]}=v{4}; end
        clear v
    end
end




General_locf=temp;
% qidssum=[temp(:,1) temp(:,contains(temp.Properties.VariableNames,{'QIDS_sum'}))];
% QIDSSuIdea=[temp(:,1) temp(:,contains(temp.Properties.VariableNames,{'QIDSSuIdea'}))];
% temp2=[QIDSSuIdea qidssum(:,2:end)];
% QIDS_locf=temp2;
writetable(General_locf,'D:\WORKINGset-D\STAR_d_suicidality\General_locf.xlsx');

clin1=temp;
%% Quality assurance - counting number of drug occurances
% m=temp; m=temp_bak
sum(sum(strcmpi(m{:,4:7},'Citalopram'))) %18030
sum(sum(strcmpi(m{:,4:7},'lithium'))) %254
sum(sum(strcmpi(m{:,4:7},'Bupropion'))) %2053
sum(sum(strcmpi(m{:,4:7},'Venlafaxine'))) %1305
sum(sum(strcmpi(m{:,4:7},'Buspirone'))) %1030
sum(sum(strcmpi(m{:,4:7},'Sertraline'))) %936


J=temp_bak;
J=J(:,~varfun(@isnumeric,temp,'output','uniform'));
sum(sum(strcmpi(J{:,:},'Citalopram'))) %17987
sum(sum(strcmpi(J{:,:},'Lithium'))) %251
sum(sum(strcmpi(J{:,:},'Bupropion')))%1964
sum(sum(strcmpi(J{:,:},'Venlafaxine')))%1237
sum(sum(strcmpi(J{:,:},'Buspirone'))) %1006
sum(sum(strcmpi(J{:,:},'Sertraline')))%879

%% QIDS LOCF
% clear all

opts = detectImportOptions('D:\WORKINGset-D\STAR_d_suicidality\QIDS2.xlsx');
varchange=opts.VariableNames([contains(opts.VariableNames,{'medication'}) & contains(opts.VariableNames,{'name'})]);
varchange2=opts.VariableNames([contains(opts.VariableNames,{'QIDS_sum'}) | contains(opts.VariableNames,{'QIDSSuIdea'}) | contains(opts.VariableNames,{'dosage'})]);

opts = setvartype(opts,varchange,'string');
opts = setvartype(opts,varchange2,'double');

QIDS=readtable('D:\WORKINGset-D\STAR_d_suicidality\QIDS2.xlsx',opts);
Qids_bak=QIDS;

level={}; level(1,1)= {'Level1_' }; level(2,1)={'Level2_'} ; level(3,1)={'Level2A_'} ; level(4,1)={'Level3_'} ; level(5,1)={'Level4_'} ; level(6,1)={'zFollow_Up_'};
temp=QIDS;

levelind(1,:)=contains(temp.Properties.VariableNames,level(1,1));
levelind(2,:)=contains(temp.Properties.VariableNames,level(2,1));
levelind(3,:)=contains(temp.Properties.VariableNames,level(3,1));
levelind(4,:)=contains(temp.Properties.VariableNames,level(4,1));
levelind(5,:)=contains(temp.Properties.VariableNames,level(5,1));
levelind(6,:)=contains(temp.Properties.VariableNames,level(6,1));

qidssum=contains(temp.Properties.VariableNames,{'QIDS_sum'});
QIDSSuIdea=contains(temp.Properties.VariableNames,{'QIDSSuIdea'});
clin=contains(temp.Properties.VariableNames,{'_Clin'});
self=contains(temp.Properties.VariableNames,{'_Self'});
med_dose1= [contains(temp.Properties.VariableNames,{'medication1'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose2= [contains(temp.Properties.VariableNames,{'medication2'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose3= [contains(temp.Properties.VariableNames,{'medication3'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose4= [contains(temp.Properties.VariableNames,{'medication4'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_name1= [contains(temp.Properties.VariableNames,{'medication1'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name2= [contains(temp.Properties.VariableNames,{'medication2'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name3= [contains(temp.Properties.VariableNames,{'medication3'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name4= [contains(temp.Properties.VariableNames,{'medication4'}) & contains(temp.Properties.VariableNames,{'name'})];



for i=1:size(temp,1)
    
    for k=1:size(levelind,1)
        
        clear locq
        locq{1,:}=temp{i,[levelind(k,:) & qidssum & clin]};
        %temp(i,[levelind(k,:) & qidssum & clin]);
        locq{2,:}=temp{i,[levelind(k,:) & qidssum & self]};
        %temp(i,[levelind(k,:) & qidssum & self]);
        locq{3,:}=temp{i,[levelind(k,:) & QIDSSuIdea & clin]};
        %temp(i,[levelind(k,:) & QIDSSuIdea & clin]);
        locq{4,:}=temp{i,[levelind(k,:) & QIDSSuIdea & self]};
        %temp(i,[levelind(k,:) & QIDSSuIdea & self]);
        v=LOCF(locq);
        clear locq
        if ~isempty(v{1}), temp{i,[levelind(k,:) & qidssum & clin]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & qidssum & self]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & QIDSSuIdea & clin]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & QIDSSuIdea & self]}=v{4}; end
        clear v
        %
        clear locd
        locd{1,:}=(temp{i,[levelind(k,:) & med_dose1]});
        locd{2,:}=(temp{i,[levelind(k,:) & med_dose2]});
        locd{3,:}=(temp{i,[levelind(k,:) & med_dose3]});
        locd{4,:}=(temp{i,[levelind(k,:) & med_dose4]});
        v=LOCF(locd);
        if ~isempty(v{1}), temp{i,[levelind(k,:) & med_dose1]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & med_dose2]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & med_dose3]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & med_dose4]}=v{4}; end
        clear v
        
        %
        clear locn
        locn{1,:}=temp{i,[levelind(k,:) & med_name1]};%table2cell(temp(i,[levelind(k,:) & med_name1]));
        locn{2,:}=temp{i,[levelind(k,:) & med_name2]};%table2cell(temp(i,[levelind(k,:) & med_name2]));
        locn{3,:}=temp{i,[levelind(k,:) & med_name3]};%table2cell(temp(i,[levelind(k,:) & med_name3]));
        locn{4,:}=temp{i,[levelind(k,:) & med_name4]};%table2cell(temp(i,[levelind(k,:) & med_name4]));
        v=LOCF(locn);
        
        
        if ~isempty(v{1}), temp{i,[levelind(k,:) & med_name1]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & med_name2]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & med_name3]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & med_name4]}=v{4}; end
        clear v
    end
end




QIDS_locf=temp;
% qidssum=[temp(:,1) temp(:,contains(temp.Properties.VariableNames,{'QIDS_sum'}))];
% QIDSSuIdea=[temp(:,1) temp(:,contains(temp.Properties.VariableNames,{'QIDSSuIdea'}))];
% temp2=[QIDSSuIdea qidssum(:,2:end)];
% QIDS_locf=temp2;
writetable(QIDS_locf,'D:\WORKINGset-D\STAR_d_suicidality\QIDS_locf3.xlsx');

clin1=temp;
%% QIDS count NANs in each level
QIDS=readtable('D:\WORKINGset-D\STAR_d_suicidality\QIDS.xlsx');

level={}; level(1,1)= {'Level1' }; level(2,1)={'Level2_'} ; level(3,1)={'Level2A'} ; level(4,1)={'Level3'} ; level(5,1)={'Level4'} ; level(6,1)={'zFollow_Up'};
temp=QIDS;
qidssum=contains(temp.Properties.VariableNames,{'QIDS_sum'});
QIDSSuIdea=contains(temp.Properties.VariableNames,{'QIDSSuIdea'});
QQ=[qidssum ; QIDSSuIdea];  QQnames={'qidssum' ; 'QIDSSuIdea'};
clsl=[contains(temp.Properties.VariableNames,{'Clin'}); contains(temp.Properties.VariableNames,{'Self'})];
clslnames={'clin' ; 'self'};
levelind(1,:)=contains(temp.Properties.VariableNames,level(1,1));
levelind(2,:)=contains(temp.Properties.VariableNames,level(2,1));
levelind(3,:)=contains(temp.Properties.VariableNames,level(3,1));
levelind(4,:)=contains(temp.Properties.VariableNames,level(4,1));
levelind(5,:)=contains(temp.Properties.VariableNames,level(5,1));
levelind(6,:)=contains(temp.Properties.VariableNames,level(6,1));
nancount=table;
for j=1:size(QQ,1)
    for k=1:size(clsl,1)
        for m=1:size(levelind,1)
            totsize=(size(temp{:,QQ(j,:) & clsl(k,:) & levelind(m,:)},1).*size(temp{:,QQ(j,:) & clsl(k,:) & levelind(m,:)},2));
            nancount(:,size(nancount,2)+1)={sum(sum(isnan(temp{:,QQ(j,:) & clsl(k,:) & levelind(m,:)})))./totsize};
            nancount.Properties.VariableNames(1,size(nancount,2))={strrep(strjoin([QQnames(j,:) clslnames(k,:) level(m,:)]),' ','_')};
            nancount(2,size(nancount,2))={totsize};
            nancount(3,size(nancount,2))={sum(sum(isnan(temp{:,QQ(j,:) & clsl(k,:) & levelind(m,:)})))};
        end
    end
end
writetable(nancount,'D:\WORKINGset-D\STAR_d_suicidality\QIDS_nancount.xlsx');
%% count nans per coloumn
HRSD=readtable('D:\WORKINGset-D\STAR_d_suicidality\HRSD.xlsx');
temp=HRSD;%QIDS;
nancnt=temp(1:2,:);
ll=size(QIDS,1);
for i=2:size(temp,2)
    nancnt(1,i)={sum(isnan(QIDS{:,i}))};
    nancnt(2,i)={sum(isnan(QIDS{:,i}))/ll};
end
nancnt(:,1)=[];
writetable(nancnt,'D:\WORKINGset-D\STAR_d_suicidality\HRSD_nancount2.xlsx');
%% MED
c=readtable('D:\WORKINGset-D\STAR_d_suicidality\3_MEDHX_ccv011.xlsx');
c.Properties.VariableNames{'src_subject_id'} = 'SubjectID';
c.Properties.VariableNames{'level'} = 'TreatmentLevel';
q=readtable('D:\WORKINGset-D\STAR_d_suicidality\QIDS_locf.xlsx');
%c=table2dataset(c);
temp=QIDS(:,[1 3 4]);
%ccc=unstack(c,[4:35],1,'AggregationFunction',@mean)
%warning off;
for i=[unique(c.SubjectID)]'
    Treat= unique(c(c.SubjectID==i,:).TreatmentLevel);
    for k=[1:size(unique(c(c.SubjectID==i,:).TreatmentLevel),1)]
        
        Day= unique(c([c.SubjectID==i & strcmpi(c.TreatmentLevel,Treat(k))],:).days_baseline);
        for j=[1:size(unique(c([c.SubjectID==i & strcmpi(c.TreatmentLevel,Treat(k))],:).days_baseline),1)]
            
            %for l=[1:size(unique(b([b.SubjectID==i & strcmpi(b.TreatmentLevel,Treat(k)) & b.DaysSinceBaseline==day(j) ],:).FormUsed_assessmentName),1)]
            %form=unique(b([b.SubjectID==i & strcmpi(b.TreatmentLevel,Treat(k)) & b.DaysSinceBaseline==day(j) ],:).FormUsed_assessmentName);
            %piv=b([b.SubjectID==i & strcmpi(b.TreatmentLevel,Treat(k)) & b.DaysSinceBaseline==day(j) & strcmpi(b.FormUsed_assessmentName,form(l))],[1 3 6 7] );
            piv=c([c.SubjectID==i & strcmpi(c.TreatmentLevel,Treat(k)) & c.days_baseline==Day(j)],[1 2 4:35] );
            
            %dou=contains(varfun(@class,piv,'OutputFormat','cell'),{'double'});
            %piv(:,dou)=cellstr(string(piv{:,dou}))
            %table2cell
            %num2str(piv(
            %piv=c([c.SubjectID==i & strcmpi(c.TreatmentLevel,Treat(k)) & c.days_baseline==day(j)],[1 2 5 9 35] );
            
            %ll=1;
            if size(piv,1)>1
                ll=[((min(sum(isnan(piv{:,[8:15, 18:34]})')'))==(sum(isnan(piv{:,[8:15, 18:34]})')'))==0];
                %ll=[((min(sum(strcmp(piv{:,[4:end]},'')')'))==(sum(strcmp(piv{:,[4:end]},'')')'))==0];
                %ll=[((min(sum(strcmp(piv{:,[4:end]},'')')'))==(sum(strcmp(piv{:,[4:end]},'')')'))==0];
                if sum(ll)<1
                    ll(2:end,:)=1;
                end
                piv(ll,:)=[];
            end
            
            if sum(temp.SubjectID==i)<1
                temp=[num2cell([i nan(1,size(temp,2)-1)]) ; temp];
            end
            
            %var1= strrep(strrep(['QIDSSuIdea_' char(Treat(k)) '_' num2str(j) '_' form{[1]}(1:4)],' ','_'),'-','_');
            %var2= strrep(strrep(['QIDS_sum_' char(Treat(k)) '_' num2str(j) '_' form{[1]}(1:4)],' ','_'),'-','_');
            %var1= strrep(strrep(['FirstMed_' char(Treat(k)) '_wk' num2str(j) '_' ],' ','_'),'-','_');
            var1= strrep([char(Treat(k)) '_day' num2str(j)],' ','_');
            piv.Properties.VariableNames(:,3:end)=cellstr(strrep(join([string([repmat([var1],size(piv(:,3:end),2),1 )]), piv.Properties.VariableNames(3:end)'])',' ','_'));
            
            %             vari{1}= strrep(strrep(['PrescribedMedCode1_' char(Treat(k)) '_wk' num2str(j) '_' ],' ','_'),'-','_');
            %             vari{2}= strrep(strrep(['PrescribedMedCode2_' char(Treat(k)) '_wk' num2str(j) '_' ],' ','_'),'-','_');
            %             vari{3}= strrep(strrep(['PrescribedMedCode3_' char(Treat(k)) '_wk' num2str(j) '_' ],' ','_'),'-','_');
            %             vari{4}= strrep(strrep(['PrescribedMedCode4_' char(Treat(k)) '_wk' num2str(j) '_' ],' ','_'),'-','_');
            %             vari{5}= strrep(strrep(['PrescribedMedCode5_' char(Treat(k)) '_wk' num2str(j) '_' ],' ','_'),'-','_');
            %             vari{6}= strrep(strrep(['PrescribedMedCode6_' char(Treat(k)) '_wk' num2str(j) '_' ],' ','_'),'-','_');
            
            if sum(contains(temp.Properties.VariableNames,[var1 '_']))==0
                %    temp=QIDS(:,[1 3 4]);
                %   temp(temp.SubjectID==i,4:)= piv(:,3:end)
                %temp1=cell(size(temp,1),size(piv(1,3:end),2));
                %for 1:size(piv,2)-2
                temp1=cell(size(temp,1),size(piv(:,3:end),2));
                temp1(temp.SubjectID==i,:)=table2cell(piv(:,3:end));
                %temp1(temp.SubjectID==i,1)= {piv{:,:}}
                temp=[temp temp1];
                %temp1=piv(1,3:end); temp1=[temp1 ; cell(size(temp,1)-1,size(piv(1,3:end),2))];
                %temp=[temp; cell(size(temp,1),size(piv,2))
                %temp1.Properties.VarNames=piv.Properties.VariableNames;
                %temp(temp.SubjectID==i,contains(temp.Properties.VariableNames,var1))= piv(:,p)
                %temp(temp.SubjectID==i,contains(temp.Properties.VariableNames,var1))= piv(:,p)
                %temp=[temp temp1];
                temp.Properties.VariableNames(:,contains(temp.Properties.VariableNames,'Var'))=(piv(:,3:end).Properties.VariableNames');
                %temp(size(temp,2)+1:size(temp,2)+1+size(piv,2))=nan(size(temp,1),size(piv,2))
                %temp.temp1=cell(size(temp,1),1);
                %temp([temp.SubjectID==i],:).temp1=piv.NameFirstMedicationThatTheParticipantHasTaken;
                %temp([temp.SubjectID==i],:)= [ temp([temp.SubjectID==i],:) piv(:,[3:end])]
                
                %temp.Properties.VariableNames{'temp1'} = var1;
            else
                %contains(varfun(@class,piv,'OutputFormat','cell'),{'cell'})
                %temp([temp.SubjectID==i],contains(temp.Properties.VariableNames,var1))=cell2table(piv(:,3:end))
                %temp([temp.SubjectID==i],contains(temp.Properties.VariableNames,var1))=table2cell(piv(:,3:end))
                temp1=table2cell(temp(:,contains(temp.Properties.VariableNames,[var1 '_'])));
                temp1([temp.SubjectID==i],:)=table2cell(piv(:,3:end));
                temp(:,contains(temp.Properties.VariableNames,[var1 '_']))=temp1;
                %table2cell(piv(:,3:end))
                %outerjoin(temp([temp.SubjectID==i],logical([1 contains(temp.Properties.VariableNames(2:end),var1)])),piv(:,[1 3:end]),'MergeKeys',true);
                %Final2=outerjoin(Final2,stats,'MergeKeys',true);
                
            end
            %             for g=1:size(vari,2)
            %                 if sum(strcmpi(temp.Properties.VariableNames,vari{g}))==0
            %                     temp.temp2=cell(size(temp,1),1);
            %                     temp([temp.SubjectID==i],:).temp2=piv{:,strcmp(piv.Properties.VariableNames,vari{g}(1:18))};
            %                     temp.Properties.VariableNames{'temp2'} = vari{g};
            %                 else
            %                     temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,vari{g}(1:18)))=piv{:,strcmp(piv.Properties.VariableNames,vari{g}(1:18))};
            %                 end
            %             end
            
            %clear -vars form temp
            %end
            %clear -vars day
        end
        %clear -vars Treat
    end
end

warning on;
bak=temp;
temp=sortrows(temp,1);
temp=[temp(:,[1:3]) temp(:,sort_nat(temp.Properties.VariableNames(4:end)))];
clin=temp(:,logical([ones(1,11) contains(temp.Properties.VariableNames(12:end),'_Clin')]));
self=temp(:,logical([1 zeros(1,10) contains(temp.Properties.VariableNames(12:end),'_Self')]));
QIDS=outerjoin(clin,self,'MergeKeys',true);
writetable(temp,'MEDtemp.xlsx')
%% HRSD

b=readtable('D:\WORKINGset-D\STAR_d_suicidality\hrsd01.xlsx'); b_bak=b;
% b=b_bak;
b=b(:,[5 6 7 11:30]);
%unique(b.TreatmentLevel)
b.Properties.VariableNames{23} = 'TreatmentLevel'; b.Properties.VariableNames{1}='SubjectID';
b(strcmpi(b.TreatmentLevel, 'Level 1'),:).TreatmentLevel=repmat({'Level1'},size(b(strcmpi(b.TreatmentLevel, 'Level 1'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 2A'),:).TreatmentLevel=repmat({'Level2A'},size(b(strcmpi(b.TreatmentLevel, 'Level 2A'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 2'),:).TreatmentLevel=repmat({'Level2'},size(b(strcmpi(b.TreatmentLevel, 'Level 2'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 3'),:).TreatmentLevel=repmat({'Level3'},size(b(strcmpi(b.TreatmentLevel, 'Level 3'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Level 4'),:).TreatmentLevel=repmat({'Level4'},size(b(strcmpi(b.TreatmentLevel, 'Level 4'),:),1),1);
b(strcmpi(b.TreatmentLevel, 'Follow Up'),:).TreatmentLevel=repmat({'zFollow-Up'},size(b(strcmpi(b.TreatmentLevel, 'Follow Up'),:),1),1);
b(strcmp(b.TreatmentLevel,''),:)=[];
b=sortrows(b,22); b=sortrows(b,23); b=sortrows(b,1);

warning off;
for i=1:size(b,1)%[unique(b.SubjectID)]' %i=2846
    %if isnan(b.QIDSTotalScore(i))
    b.HRSTotalScore_recorded_(i)=sum(b{i,4:20});
    %max(b{i,6:9})+max(b{i,11:14})+max(b{i,20:21})+b{i,10}+nansum(b{i,15:19});
    %end
end
d=b;
%% HRSD
%writetable(d,'D:\WORKINGset-D\STAR_d_suicidality\hrsd02.xlsx');
d=readtable('D:\WORKINGset-D\STAR_d_suicidality\hrsd02.xlsx');
m=readtable('D:\WORKINGset-D\STAR_d_suicidality\STARD-medications.csv','Format','%f %s %f %s %s %s %s %f %f %f %f');
d_backup=d;  m_backup=m;
%d=d(:,[1 2 3 15 21 22 23]);
%temp=QIDS(:,[1 3 4]);
temp=table; temp.SubjectID=unique(d.SubjectID);
%warning off;
for i=[unique(d.SubjectID)]'
    Treat= unique(d(d.SubjectID==i,:).TreatmentLevel);
    for k=[1:size(unique(d(d.SubjectID==i,:).TreatmentLevel),1)]
           Day= sort(unique(d([d.SubjectID==i & strcmpi(d.TreatmentLevel,Treat(k))],:).DaysSinceBaseline));
           %days1=[1:size(Day,1)];
           days1=[1:size(unique(d([d.SubjectID==i & strcmpi(d.TreatmentLevel,Treat(k))],:).DaysSinceBaseline),1)];
        for j=days1
            piv=table;
            piv=d([d.SubjectID==i & strcmpi(d.TreatmentLevel,Treat(k)) & d.DaysSinceBaseline==Day(j)],[1 3 4:22] );
            
            if size(piv,1)>1
                ll=[((min(sum(isnan(piv{:,[3 4]})')'))==(sum(isnan(piv{:,[3 4]})')'))==0];
                if sum(ll)<1
                    ll(1:end-1,:)=1;
                end
                piv(ll,:)=[];
                if size(piv,1)>1; piv=piv(end,:); end
            end
            if sum(temp.SubjectID==i)<1
                temp=[num2cell([i nan(1,size(temp,2)-1)]) ; temp];
            end
            %Week=sort(round(Day/7));
            %weeknum=1:size(Week,1);
            Week=sort(floor(Day/7)-floor(Day(1)/7));
                if Week(1)==0
                    weeknum=[0:size(Week,1)-1];
                else
                    weeknum=[1:size(Week,1)]; % weeknum=1:size(Week,1);
                end
                
            var1= strrep(strrep(['HRSSuicide_' char(Treat(k)) '_week' num2str(weeknum(j)) '_'  ],' ','_'),'-','_');
            var2= strrep(strrep(['HRSTotalScore_' char(Treat(k)) '_week' num2str(weeknum(j)) '_' ],' ','_'),'-','_');
            
            if sum(strcmpi(temp.Properties.VariableNames,var1))==0
                temp.temp1=nan(size(temp,1),1);
                temp([temp.SubjectID==i],:).temp1=cell2mat({piv.HRSSuicide});
                temp.Properties.VariableNames{'temp1'} = var1;
            else
                temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,var1))={piv.HRSSuicide};
            end
            if sum(strcmpi(temp.Properties.VariableNames,var2))==0
                temp.temp2=nan(size(temp,1),1);
                temp([temp.SubjectID==i],:).temp2=cell2mat({piv.HRSTotalScore_recorded_});
                temp.Properties.VariableNames{'temp2'} = var2;
            else
                temp([temp.SubjectID==i],strcmpi(temp.Properties.VariableNames,var2))={piv.HRSTotalScore_recorded_};
            end
            %clear piv
            pivm=table;
            med_wk=sort(m{[m.SubjectID==i & strcmpi(m.TreatmentLevel,Treat(k)) ],3},'ascend');
            if ~isempty(med_wk)
            mesure_dist=pdist2(med_wk,Week);
            
            if j==days1(end)
                gap1=4.1;
            elseif j==days1(1) && Week(j)==0 && med_wk(1)==0
                gap1=0.1;
            else
                gap1=1.1;
            end
            if min(mesure_dist(:,j))<gap1
            [yo, II]=min(mesure_dist(:,j));
            pivm=sortrows([ m([m.SubjectID==i & strcmpi(m.TreatmentLevel,Treat(k)) ],[3:end])],1);
            % pivm=m([m.SubjectID==i & strcmpi(m.TreatmentLevel,Treat(k)) & m.Med_week==Week(j)],[3:end] );
            pivm=pivm(II,:);
            end
            end
            
            if size(pivm,1)>1
                [~,loc5]=min(sum(cellfun(@isempty,(table2cell(pivm(:,[1 5])))),2));
                pivm=pivm(loc5,:);
                if size(pivm,1)>1; pivm=pivm(1,:); end
            end
            if ~isempty(pivm)
                pivm=[piv(:,21) pivm];
%             else
%                 pivm=piv(:,21);
%             end
            %if ~isempty(pivm)
                var3=['_' char(Treat(k)) '_week' num2str(weeknum(j)) '_'];
                pivm.Properties.VariableNames=strrep(strcat(pivm.Properties.VariableNames,var3),'-','_');
                if sum(contains(temp.Properties.VariableNames,pivm.Properties.VariableNames(1)))==0
                    temp3=num2cell(nan(size(temp,1),size(pivm,2))); %array2table(nan(size(temp,1),size(pivm,2)));
                    temp3([temp.SubjectID==i],:)=table2cell(pivm);
                    temp3=array2table(temp3);
                    temp3.Properties.VariableNames=pivm.Properties.VariableNames;
                    temp=[temp temp3];
                    clear temp3
                else
                    temp{[temp.SubjectID==i],pivm.Properties.VariableNames}=table2cell(pivm);
                end
            end
               % clear pivm piv
            end
        end
        
    end
end
bak=temp;
writetable(temp,'D:\WORKINGset-D\STAR_d_suicidality\HRSDtemp2.xlsx')

%warning on;
temp=sortrows(temp,1);
wk0=contains(temp.Properties.VariableNames(1:end),'week0_');
wk1=contains(temp.Properties.VariableNames(1:end),'week1_');
wk2=contains(temp.Properties.VariableNames(1:end),'week2_');
wk3=contains(temp.Properties.VariableNames(1:end),'week3_');
wk4=contains(temp.Properties.VariableNames(1:end),'week4_');
wk5=contains(temp.Properties.VariableNames(1:end),'week5_');

temp=[temp(:,[1]) temp(:,[wk0]) temp(:,[wk1]) temp(:,[wk2]) temp(:,[wk3]) temp(:,[wk4]) temp(:,[wk5])];



%temp=[temp(:,[1:3]) temp(:,sort_nat(temp.Properties.VariableNames(4:end)))];
lvl0=contains(temp.Properties.VariableNames(1:end),'Enrollment');
lvl1=contains(temp.Properties.VariableNames(1:end),'Level1_');
lvl2=contains(temp.Properties.VariableNames(1:end),'Level2_');
lvl2A=contains(temp.Properties.VariableNames(1:end),'Level2A');
lvl3=contains(temp.Properties.VariableNames(1:end),'Level3');
lvl4=contains(temp.Properties.VariableNames(1:end),'Level4');
flup=contains(temp.Properties.VariableNames(1:end),'zFollow_Up');
temp=[temp(:,[1]) temp(:,[lvl0]) temp(:,[lvl1]) temp(:,[lvl2]) temp(:,[lvl2A]) temp(:,[lvl3]) temp(:,[lvl4]) temp(:,[flup])];

Suicide=temp(:,logical([1 contains(temp.Properties.VariableNames(2:end),'HRSSuicide_')]));
Total=temp(:,logical([1 contains(temp.Properties.VariableNames(2:end),'HRSTotalScore_')]));
HRSD=outerjoin(Suicide,Total,'MergeKeys',true);
writetable(temp,'D:\WORKINGset-D\STAR_d_suicidality\HRSDtemp.xlsx')
%writetable(HRSD,'D:\WORKINGset-D\STAR_d_suicidality\HRSD.xlsx')
%% HRSD_locf
opts = detectImportOptions('D:\WORKINGset-D\STAR_d_suicidality\HRSDtemp.xlsx');
varchange=opts.VariableNames([contains(opts.VariableNames,{'medication'}) & contains(opts.VariableNames,{'name'})]);
varchange2=opts.VariableNames([contains(opts.VariableNames,{'HRSSuicide'}) | contains(opts.VariableNames,{'HRSTotalScore'}) | contains(opts.VariableNames,{'dosage'})]);
%[varchange' varchange2']
opts = setvartype(opts,varchange,'string');
opts = setvartype(opts,varchange2,'double');
HRSD=readtable('D:\WORKINGset-D\STAR_d_suicidality\HRSDtemp.xlsx',opts);
HRSD_bak=HRSD;
temp=HRSD;


level={}; level(1,1)= {'Enrollment' }; level(2,1)= {'Level1_' }; level(3,1)={'Level2_'} ; level(4,1)={'Level2A_'} ; level(5,1)={'Level3_'} ; level(6,1)={'Level4_'} ; level(7,1)={'zFollow_Up_'};
levelind(1,:)=contains(temp.Properties.VariableNames,level(1,1));
levelind(2,:)=contains(temp.Properties.VariableNames,level(2,1));
levelind(3,:)=contains(temp.Properties.VariableNames,level(3,1));
levelind(4,:)=contains(temp.Properties.VariableNames,level(4,1));
levelind(5,:)=contains(temp.Properties.VariableNames,level(5,1));
levelind(6,:)=contains(temp.Properties.VariableNames,level(6,1));
levelind(7,:)=contains(temp.Properties.VariableNames,level(7,1));

HRSsui=contains(temp.Properties.VariableNames,{'HRSSuicide'});
HRStot=contains(temp.Properties.VariableNames,{'HRSTotalScore'});

med_dose1= [contains(temp.Properties.VariableNames,{'medication1'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose2= [contains(temp.Properties.VariableNames,{'medication2'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose3= [contains(temp.Properties.VariableNames,{'medication3'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_dose4= [contains(temp.Properties.VariableNames,{'medication4'}) & contains(temp.Properties.VariableNames,{'dosage'})];
med_name1= [contains(temp.Properties.VariableNames,{'medication1'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name2= [contains(temp.Properties.VariableNames,{'medication2'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name3= [contains(temp.Properties.VariableNames,{'medication3'}) & contains(temp.Properties.VariableNames,{'name'})];
med_name4= [contains(temp.Properties.VariableNames,{'medication4'}) & contains(temp.Properties.VariableNames,{'name'})];

for i=1:size(temp,1)
    for k=1:size(levelind,1)
        clear locq
        locq{1,:}=temp{i,[levelind(k,:) & HRSsui]};
        locq{2,:}=temp{i,[levelind(k,:) & HRStot]};
        v=LOCF(locq);
        clear locq
        if ~isempty(v{1}), temp{i,[levelind(k,:) & HRSsui]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & HRStot]}=v{2}; end
        clear v
        clear locd
        locd{1,:}=(temp{i,[levelind(k,:) & med_dose1]});
        locd{2,:}=(temp{i,[levelind(k,:) & med_dose2]});
        locd{3,:}=(temp{i,[levelind(k,:) & med_dose3]});
        locd{4,:}=(temp{i,[levelind(k,:) & med_dose4]});
        v=LOCF(locd);
        if ~isempty(v{1}), temp{i,[levelind(k,:) & med_dose1]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & med_dose2]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & med_dose3]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & med_dose4]}=v{4}; end
        clear v
        %
        clear locn
        locn{1,:}=temp{i,[levelind(k,:) & med_name1]};
        locn{2,:}=temp{i,[levelind(k,:) & med_name2]};
        locn{3,:}=temp{i,[levelind(k,:) & med_name3]};
        locn{4,:}=temp{i,[levelind(k,:) & med_name4]};
        v=LOCF(locn);
        if ~isempty(v{1}), temp{i,[levelind(k,:) & med_name1]}=v{1}; end
        if ~isempty(v{2}), temp{i,[levelind(k,:) & med_name2]}=v{2}; end
        if ~isempty(v{3}), temp{i,[levelind(k,:) & med_name3]}=v{3}; end
        if ~isempty(v{4}), temp{i,[levelind(k,:) & med_name4]}=v{4}; end
        clear v
    end
end


HRSD_locf=temp;
writetable(HRSD_locf,'D:\WORKINGset-D\STAR_d_suicidality\HRSD_locf.xlsx');
% for i=1:size(temp,1)
%     for j=1:size(QQ,1)
%
%         loc1=find(diff(isnan(temp{i,[QQ(j,:)  ]})))+1;
%         loc2=loc1(1:2:end);
%         loc3=[loc1(2:2:end) size(temp{i,[QQ(j,:) ]},2)+1];
%         if ~isempty(loc1)
%             for ff=1:size(loc2,2)
%                 temp{i,[QQ(j,:) ]}(:,loc2(ff):loc3(ff)-1)=temp{i,[QQ(j,:) ]}(:,[loc2(ff)-1]);
%             end
%         end
%
%     end
% end

%% HRSD count NANs in each level
level={}; level(1,1)= {'Enrollment' }; level(2,1)= {'Level1' }; level(3,1)={'Level2_'} ; level(4,1)={'Level2A'} ; level(5,1)={'Level3'} ; level(6,1)={'Level4'} ; level(7,1)={'zFollow_Up'};
temp=HRSD;
qidssum=contains(temp.Properties.VariableNames,{'HRSTotalScore_'});
QIDSSuIdea=contains(temp.Properties.VariableNames,{'HRSSuicide_'});
QQ=[qidssum ; QIDSSuIdea];  QQnames={'HRSTotalScore_' ; 'HRSSuicide_'};
%clsl=[contains(temp.Properties.VariableNames,{'Clin'}); contains(temp.Properties.VariableNames,{'Self'})];
%clslnames={'clin' ; 'self'};
levelind(1,:)=contains(temp.Properties.VariableNames,level(1,1));
levelind(2,:)=contains(temp.Properties.VariableNames,level(2,1));
levelind(3,:)=contains(temp.Properties.VariableNames,level(3,1));
levelind(4,:)=contains(temp.Properties.VariableNames,level(4,1));
levelind(5,:)=contains(temp.Properties.VariableNames,level(5,1));
levelind(6,:)=contains(temp.Properties.VariableNames,level(6,1));
levelind(7,:)=contains(temp.Properties.VariableNames,level(7,1));
nancount=table;
for j=1:size(QQ,1)
    
    for m=1:size(levelind,1)
        totsize=(size(temp{:,QQ(j,:) & levelind(m,:)},1).*size(temp{:,QQ(j,:) & levelind(m,:)},2));
        nancount(:,size(nancount,2)+1)={sum(sum(isnan(temp{:,QQ(j,:) & levelind(m,:)})))./totsize};
        nancount.Properties.VariableNames(1,size(nancount,2))={strrep(strjoin([QQnames(j,:) level(m,:)]),' ','_')};
        nancount(2,size(nancount,2))={totsize};
        nancount(3,size(nancount,2))={sum(sum(isnan(temp{:,QQ(j,:) & levelind(m,:)})))};
    end
    
end
writetable(nancount,'D:\WORKINGset-D\STAR_d_suicidality\HRSD_nancount.xlsx');
%%
st='Enrollment'
Enrollment=a(strcmpi(a.TreatmentLevel,'Enrollment'),:);
Enrollment=outerjoin(a(strcmpi(a.TreatmentLevel,st),:),b(strcmpi(b.TreatmentLevel,st),[1 4 6 7]),'MergeKeys',true);
%Enrollment=outerjoin(Enrollment,c(strcmpi(c.TreatmentLevel,st),[1 4 5 6 7 8 9 10 11 12]),'MergeKeys',true);
Enrollment=outerjoin(Enrollment,d(strcmpi(d.TreatmentLevel,st),[1 4 5]),'MergeKeys',true);
%%
st='Follow up'



fff=a(strcmpi(a.TreatmentLevel,st),:);
fff=sortrows(fff,2);
fff=sortrows(fff,1);



%base
%base
% (ismember(unique(temp.SubjectID),unique(b.SubjectID)))
% size(unique(b.SubjectID))
% size(unique(temp.SubjectID))

base.SubjectID==1
subsT.Properties.VariableNames{'ISP'} = strrep(strjoin({'ISP' strrep(strjoin([elecname elecHOMOname]),' ','_') strrep(num2str([statwin(1) statwin(2)]),'  ','_')}, '_'),'__','_');




DaysSinceBaseline
for j=1:size(fff.TreatmentLevel(fff.SubjectID==i),1)
    fff(fff.SubjectID==i,:).TreatmentLevel(j)={['Follow up_'  num2str(j)]};
end


ff1=fff([1
    unstack(fff([1,'TreatmentLevel',[2])



unique(b.TreatmentLevel)
sortrows(b,1);
st='Follow up_1'
Follow_up=outerjoin(fff(strcmpi(fff.TreatmentLevel,st),:),b(strcmpi(b.TreatmentLevel,'Follow-Up'),[1 4 6 7]),'MergeKeys',true);
Follow_up=outerjoin(Follow_up,c(strcmpi(c.TreatmentLevel,st),[1 4 5 6 7 8 9 10 11 12]),'MergeKeys',true);
Follow_up=outerjoin(Follow_up,d(strcmpi(d.TreatmentLevel,st),[1 4 5]),'MergeKeys',true);
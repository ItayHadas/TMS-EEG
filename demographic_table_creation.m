
load('D:\WORKINGSET_D\Trial_compare_SCD_SCS_SGC_hipp\3trials_stats_struct.mat')

mean(comp2.MST.clinical.AgeAtEnrollment)
nanmean(comp2.MST.clinical.x_OfAcuteTxsReceived)
nanstd(comp2.MST.clinical.x_OfAcuteTxsReceived)
sum(strcmpi(comp2.MST.clinical.Sex,'f'))

sum(comp2.MST.HAMD.LOCF_responders) ; 14/31
nanmean(comp2.MST.MOCA.Moca_TotalBaseline)
nanstd(comp2.MST.MOCA.Moca_TotalBaseline)
nanmean(comp2.MST.MOCA.Moca_TotalAcutePost_MaintPre)
sum(~isnan(comp2.MST.MOCA.Moca_TotalAcutePost_MaintPre))

% ECT
ECT=table;
ECT.age_mean=mean(comp2.ECT.Clinical.AgeAtTxStart)
ECT.age_std=nanstd(comp2.ECT.Clinical.AgeAtTxStart)
ECT.N_Females=sum(strcmpi(comp2.ECT.Clinical.Sex,'f'))
ECT.N_males=sum(strcmpi(comp2.ECT.Clinical.Sex,'m'))
ECT.avg_Number_treatments=nanmean(comp2.ECT.Clinical.x_OfTxsReceived)
ECT.std_Number_treatments=nanstd(comp2.ECT.Clinical.x_OfTxsReceived)
ECT.Responders=sum(comp2.ECT.HAMD.LOCF_responders2)
ECT.remitters=sum(comp2.ECT.HAMD.LOCF_remitters3)
ECT.HAMD_baseline_avg=nanmean(comp2.ECT.HAMD.TotalPre)
ECT.HAMD_baseline_std=nanstd(comp2.ECT.HAMD.TotalPre)
ECT.HAMD_post_avg=nanmean(comp2.ECT.HAMD.TotalPost)
ECT.HAMD_post_std=nanstd(comp2.ECT.HAMD.TotalPost)
writetable(ECT,'ECT_demographics.xls')


sum(comp2.MST.HAMD.LOCF_responders) ; 14/31
nanmean(comp2.MST.MOCA.Moca_TotalBaseline)
nanstd(comp2.MST.MOCA.Moca_TotalBaseline)
nanmean(comp2.MST.MOCA.Moca_TotalAcutePost_MaintPre)
sum(~isnan(comp2.MST.MOCA.Moca_TotalAcutePost_MaintPre))
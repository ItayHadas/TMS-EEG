%% creating individualized table - line per participant/reccording

cd D:\WORKINGset_D\MST

[fileNames, pathName]=Z_getSetsFileNames('set');
preprocess_stats=table;
for i=1: size(fileNames,1)
    fileName=fileNames{i};
    EEG = pop_loadset( [pathName fileName]);
    
    preprocess_stats.subject(i)=str2num(fileName(4:6));
    preprocess_stats.stage{i}=pathName(21:24);
    [~, t1] = min(abs(EEG.times-(25)));
    [~, t2] = min(abs(EEG.times-(30)));
    data=[];   data=EEG.data(:,[t1:end],:);
    preprocess_stats.rank(i)=rank(double(data(:,:)));
    preprocess_stats.epochs_removed(i)=size(EEG.urevent,2)-size(data,3);
    preprocess_stats.Channels_removed(i)=size(EEG.chanlocs,2)-size(EEG.icachansind,2);
    if isfield(EEG,'comp_removed')
    preprocess_stats.ICA_component_removed(i)=size(EEG.comp_removed,2);
    else
     preprocess_stats.ICA_component_removed(i)=nan;
    end
end


preprocess_stats_pre = preprocess_stats;

preprocess_stats_post =preprocess_stats;

stats=[preprocess_stats_pre ; preprocess_stats_post]

preprocess_stats_post = readtable('preprocess_stats_post.xlsx');

%% group descriptive statistics

var=6

[h, p, ci, stats_tt]=ttest(cell2mat(preprocess_stats_pre{:,var}),cell2mat(preprocess_stats_post{:,var}))


sta = varfun(@(x) [mean(x,'omitnan'); std(x,'omitnan')],stats,'GroupingVariables','stage' ,'InputVariables',[3:6] )

%%
mean([preprocess_stats_post.epochs_removed{:}])  
mean([preprocess_stats_pre.epochs_removed{:}]) 
std([preprocess_stats_pre.epochs_removed{:}]) 
[h, p, ci, stats]=ttest([preprocess_stats_pre.epochs_removed{:}],[preprocess_stats_post.epochs_removed{:}])
[h, p, ci, stats]=ttest2([preprocess_stats_pre.epochs_removed{:}],[preprocess_stats_post.epochs_removed{:}])
[h, p, stats]=ranksum([preprocess_stats_pre.epochs_removed{:}],[preprocess_stats_post.epochs_removed{:}])

mean([preprocess_stats_pre.Channels_removed{:}]) 
std([preprocess_stats_pre.Channels_removed{:}]) 
[h, p, ci, stats]=ttest([preprocess_stats_pre.Channels_removed{:}],[preprocess_stats_post.Channels_removed{:}])
[h, p, ci, stats]=ttest2([preprocess_stats_pre.epochs_removed{:}],[preprocess_stats_post.epochs_removed{:}])
[h, p, stats]=ranksum([preprocess_stats_pre.epochs_removed{:}],[preprocess_stats_post.epochs_removed{:}])
%%
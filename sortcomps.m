function [EEG, varsPerc] = sortcomps( EEG, time) 

%Check input for ICA weights
if isempty(EEG.icawinv)
	error('No ICA weights. Please run ICA first.');
end
%time=[0 200];
[~, T1] = min(abs(EEG.times-(time(1))));
[~, T2] = min(abs(EEG.times-(time(2))));
vars = arrayfun(@(x)var(mean(eeg_getdatact(EEG, 'component', [x], 'projchan', [],'samples',[T1:T2]),3)),1:size(EEG.icawinv,2)); %extracts time course for each component
vars_norm = vars/sum(vars)*100; %calculates the % variance of each component relative to all components
[xSorted, ixsSort] = sort(vars_norm, 'descend'); %ranks components based on %variance
varsPerc = vars_norm(ixsSort); 

EEG.icawinv = EEG.icawinv(:,ixsSort); %alters inverse weights matrix based on component variance order
EEG.icaweights = EEG.icaweights(ixsSort,:); %alters weights matrix based on component variance order
if ~isempty(EEG.icaact)
    EEG.icaact = EEG.icaact(ixsSort,:,:); %alters time course matrix based on component variance order
end

fprintf('ICA weights sorted by time course variance\n');
        
end

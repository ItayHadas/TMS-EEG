function EEG = COMPremove(EEG,fn,comp)
setname=EEG.setname;
filename=[EEG.setname '.cnt'];
datfile=[EEG.setname '.fdt'];

EEG = pop_subcomp( EEG, comp, 0);


EEG.setname=[setname  fn];
EEG.filename=[EEG.setname '.cnt'];
EEG.datfile=[EEG.setname '.fdt'];
EEG.Removed_ICA_Comps.Removed_Comps2=num2cell(comp);

EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, [EEG.filename]);
end

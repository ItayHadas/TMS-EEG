% trimming and interpulating around the event in epoched EEG data.

function  EEG = Triminterp(EEG,trim)
%trim=[-6 16];
srate=EEG.srate;
trimt=trim*srate/1000;
[~, loc]=min(abs(EEG.times));
[~, basl] = min(abs(EEG.times-(-700)));
data=double(EEG.data);

for elec = 1:size(EEG.data,1)
    for epo = 1:size(EEG.data,3)
        % baseline is not reccomended in this time - strong slow frequency might create issues.
        % 
        %data(elec,[1:loc],epo)=data(elec,[1:loc],epo)-mean(data(elec,[loc-[800*srate/1000]:loc-[400*srate/1000]],epo)); %baselline correct for negative times
        %data(elec,[loc:end],epo)=data(elec,[loc:end],epo)-mean(data(elec,[loc+[600*srate/1000]:loc+[950*srate/1000]],epo)); %baselline correct for positivve times
        %data(elec,:,epo)=data(elec,:,epo)-mean(data(elec,[basl-200:basl],epo),2); %baselline correct 
        % plot(data(elec,[1:loc],epo)-mean(data(elec,[loc-[800*srate/1000]:loc-[400*srate/1000]],epo)))
        % EEG.times(find(isnan(data(elec,:,epo))))
        xx=data(elec,[loc+round(trimt(1))-25:loc+round(trimt(2))+25],epo);
        xx(25:end-26)=nan;
        inn = ~isnan(xx);
        i1 = (1:numel(xx)).';
        pp = interp1(i1(inn),xx(inn),'pchip','pp'); % spline 'linear'
        out_1 = fnval(pp,linspace(i1(1),i1(end),size(xx,2)));
        data(elec,[loc+round(trimt(1))-25:loc+round(trimt(2))+25],epo)=out_1;
        
        %figure; hold on; plot(EEG.times([loc+trimt(1)-25:loc+trimt(2)+25]),data(elec,[loc+trimt(1)-25:loc+trimt(2)+25],epo))
    end
end

EEG.data=double(data);

if sum(strcmpi(fields(EEG),'tmscut'))>0
    if   strcmpi(class(EEG.tmscut),'double')
        EEG.tmscut.cutTimesTMS=trim;
        
    elseif sum(strcmpi(fields(EEG.tmscut),'cutTimesTMS'))>0
        EEG.tmscut.cutTimesTMS2=trim;
    end
else
    EEG.tmscut.cutTimesTMS=trim;
end

EEG = eeg_checkset( EEG );
disp(['trimmed ' trim])

end
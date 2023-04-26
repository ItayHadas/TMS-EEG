function EEGdata = FilterTEP(EEGdata, low, high,srate)

%% eeglab filter ,elecs,srate
%  low=0.2;
%  high=80; EEGdata=EEG.data; elecs=[1:size(EEG.data,1)]; srate=EEG.srate;
%profile on
% filtorder=[];
% usefft=0;
% plotfreqz=0;
% firtype='firls';
% causal=0;
%revfilter=0;
%EEG = pop_eegfiltnew(EEG, low, high,  filtorder, revfilter, usefft, plotfreqz);
%addname=['_bpf' num2str(low) '_' num2str(high)];



% ch=size(elecs,2);
% %padding against filter edge effect (reza style)
% if ep==1
%     K2=[];  K2=[EEG.data(:,1:500) EEG.data  EEG.data(:,end-500:end)];
% elseif ep>1
%     K2=[];  K2=[EEG.data(:,1:500,:) EEG.data  EEG.data(:,end-500:end,:)];
% end
K2=double([EEGdata]); 
%[x_filt,y_filt]=butter(2,[low high]/(1000/2),'bandpass');

[x_filthigh,y_filthigh]=butter(4,low./(srate/2),'high'); %Highpass
%[zh, ph, kh]=butter(10,low./(srate/2),'high'); %Highpass
%[x_filthigh,y_filthigh] = zp2sos(zh,ph,kh);
[x_filtlow,y_filtlow]=butter(4,high./(srate/2)); %Lowpass
%[zl, pl, kl]=butter(10,high./(srate/2)); %Lowpass
%[x_filtlow,y_filtlow] = zp2sos(zl,pl,kl);

%  [B,A] = BUTTER(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
[x_filtstop,y_filtstop]=butter(4,[58./(srate/2) 62./(srate/2)],'stop');
%[zst, pst, kst]=butter(10,[58./(srate/2) 62./(srate/2)],'stop'); %Lowpass
%[x_filtstop,y_filtstop] = zp2sos(zst,pst,kst);
[ch,~,ep]=size(EEGdata);
if ep==1
    for elec=1:size(EEGdata,1)
        %EEG.data(elec,:)=filtfilt(x_filt,y_filt,double(EEG.data(elec,:)));
        %K2(elec,:)=filtfilt(x_filt,y_filt,double(K2(elec,:)));
        K2(elec,:)=filtfilt(x_filthigh,y_filthigh,K2(elec,:));
        K2(elec,:)=filtfilt(x_filtlow,y_filtlow,K2(elec,:));
        K2(elec,:)=filtfilt(x_filtstop,y_filtstop,K2(elec,:));
    end
elseif ep>1
    for elec=1:size(EEGdata,1)
        for epo=1:size(EEGdata,3)
            %EEG.data(elec,:)=filtfilt(x_filt,y_filt,double(EEG.data(elec,:)));
            %K2(elec,:)=filtfilt(x_filt,y_filt,double(K2(elec,:)));
            K2(elec,:,epo)=filtfilt(x_filthigh,y_filthigh,K2(elec,:,epo));
            K2(elec,:,epo)=filtfilt(x_filtlow,y_filtlow,K2(elec,:,epo));
            K2(elec,:,epo)=filtfilt(x_filtstop,y_filtstop,K2(elec,:,epo));
        end
    end
    
end


% removes the padding
% K2 = reshape(K2,ch,[],ep);
% K2=double(K2);
% if ep==1
%     EEG.data=K2(:,501:end-501);
% elseif ep>1
%     EEG.data=K2(:,501:end-501,:);
% end
EEGdata=double(K2);


% addname=['_pass' num2str(low) '_' num2str(high)];
% EEG = Z_append(EEG, addname); 
% EEG = eeg_checkset( EEG );
% % figure; pwelch(EEG.data',[] ,[], [], EEG.srate); xlim([0 0.12])
% mkdir([pathName '\FILT'])
% EEG = pop_saveset( EEG, [pathName 'FILT\' EEG.filename]);
end

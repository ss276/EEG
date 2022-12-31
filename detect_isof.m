function [stg,sup_datapoint]=detect_isof(fs,eegdata,emgdata,baseline_t,gap_sec)
% fs: sampling rate
% eegdata, emgdata : 1-row or column matrix
% baseline time (min) : time (min) before isoflurane induction 

% threshold  [thr1,thr2,thr3]
% either baseline time or threshold needed
% if threshold exist, no baseline calculating

if nargin < 4
    error('not enough data')
elseif nargin < 5
    gapsec=10; % default
end

data_size= size(eegdata);
if data_size(2) > data_size(1)
    eegdata=eegdata';
    emgdata=emgdata';
end

if length(eegdata)~=length(emgdata)
error('EEG and EMG data length need to be same')
end

data_size= size(eegdata);
data_time=data_size(1)/fs;
T=[1:length(eegdata)]/fs/60;

%% pre-processing
eegdata(:,1)=eegdata(:,1)-smooth(eegdata(:,1),fs*4);
eegdata(:,1) = zerofilt(eegdata(:,1),0.5,128,fs);
emgdata(:,1)=emgdata(:,1)-smooth(emgdata(:,1),fs*4);

%--------------------------------
% figure;plot(T,eegdata(:,1),'k')
% set(gcf,'position',[45 45 1100 200]);xlim([0 10]);ylim([-0.7 0.7])
% ylabel('mV'); xlabel('time (min)')
% figure;plot(T,emgdata(:,1),'k')
% set(gcf,'position',[45 45 1100 200]);xlim([0 10]);ylim([-8 7])
% ylabel('mV'); xlabel('time (min)')

%% Burst supression
% threshold for burst-suppression detection
bl=eegdata(1:baseline_t*60*fs,1);
thr1=mean(bl)+std(bl)*0.8
thr2=mean(bl)+std(bl)*1.1

clear bdata
bdata=abs(eegdata(:,1));
svec=zeros(1,length(bdata));
clear tmp* a b
tmp1=find(bdata<thr1);% threshold for suppression period
tmp2=find(bdata<thr2); % sub-threshold for sensitivity
svec(tmp1)=1;
a=find(diff(tmp1)>1 & diff(tmp1)<fix(fs/10)); % 0.1 s
for ind=1:length(a)
b=tmp1(a(ind)):tmp1(a(ind)+1);
if ismember(b,tmp2)
svec(b)=1;
end
end

clear tmp*
tmp1=find(svec==1);
tmp2=diff(tmp1); 
tmp3=find(tmp2>1); tmp3=[0, tmp3];
tmp4=find(diff(tmp3)>fs);% 1 second rule 
tmp5(:,1)=tmp3(tmp4)+1;
tmp5(:,2)=tmp3(tmp4+1);
tmp6(:,1)=tmp1(tmp5(:,1));
tmp6(:,2)=tmp1(tmp5(:,2));

sup_datapoint=tmp6; % suppression data point
bsmtx=zeros(length(eegdata(:,1)),1);  % suppression period vector
for ind=1:length(sup_datapoint)
bsmtx(sup_datapoint(ind,1):sup_datapoint(ind,2))=1;
end
bsmtx(1:1:baseline_t*60*fs)=0;


%--------------------------------
% figure;plot(T,eegdata(:,1),'k')
% hold on;plot(T(tmp6(:,1)),eegdata(tmp6(:,1),1),'r^','linewidth',2)
% hold on;plot(T(tmp6(:,2)),eegdata(tmp6(:,2),1),'gs','linewidth',2)
% set(gcf,'position',[45 45 1100 200]);xlim([9 9.5]);ylim([-0.3 0.3])
% figure;plot(T,bdata,'color',[1 0.8 0.8])
% yline(thr1,'k','linewidth',1.5); %grid on;
% hold on;plot(T(tmp6(:,1)),bdata(tmp6(:,1),1),'r^','linewidth',2)
% %hold on;plot(T(tmp6(:,2)),bdata(tmp6(:,2),1),'rs','linewidth',2)
% set(gcf,'position',[45 45 1100 200]);xlim([0 10]);ylim([0 0.5])
% set(gca,'xtick',[],'ytick',[])
%% spectrogram
fs=256;
clear spec_ifft;clear spec_f;clear spec_t;clear spec_psd2;clear psd; clear fft;
w_size = fix(1.5*fs);     % window size in datapoint
move_length = fix(0.15*fs);      % moving window size in datapoint
nfft = 2^nextpow2(w_size);
[spec_ifft, spec_f, spec_t, spec_psd] = spectrogram(eegdata(:,1), w_size, w_size-move_length,nfft, fs);

% figure;imagesc(spec_t,spec_f,spec_psd)
% set(gca,'YDir','normal','fontsize',16); ylim([0 20]);
% caxis([0 0.002])
% xlabel('time (s)','fontsize',18); ylabel('frequecny (Hz)','fontsize',18)
% colormap(jet)

figure
contourf(spec_t/60,spec_f,spec_psd,70,'lines','none')
%caxis([0 0.001])
set(gca,'YDir','normal','fontsize',16); ylim([0 15]);
ylabel('frequency(hz)')
set(gcf,'position',[45 45 1100 250])
%% Delta & EMG
% EMG
clear tmp* a* b3 b4 c*
emgabs=abs(emgdata);
tmp2=smooth(emgabs,10*fs);
a3=nanmedian(tmp2(1:find(T>baseline_t,1)),1);
a4=ones(length(emgabs),1); a4(find(tmp2>a3))=0;

% 2 Hz
pow1=mean(spec_psd(2:find(spec_f>=4,1),:),1); % ~2Hz
tmp3=smooth(pow1,find(spec_t>=10,1)); % ~10s
%figure;plot(spec_t,tmp3,'linewidth',2,'k')
tmp3 = interp(tmp3,move_length);
tmp3=[ones(w_size/2-1,1)*tmp3(1);tmp3;];
tmp3=[tmp3;ones(length(T)-length(tmp3),1)*tmp3(end)];
% figure;plot(T,tmp3,'k','linewidth',2)
% hold on;plot(spec_t/60,pow1)
b3=mean(tmp3(1:find(T>baseline_t,1)))+std(tmp3(1:find(T>baseline_t,1)))*2
b4=zeros(length(tmp3),1); b4(find(tmp3>b3))=1;
% hold on; plot(T,b4*0.003)

c4=a4.*b4;
c4(1:baseline_t*60*fs)=0;

%--------------------------------
% figure; plot(T,tmp2,'color',[0.8 0.8 0.5],'linewidth',2)
% yline(a3,'linewidth',1.5)
% set(gcf,'position',[45 45 1100 200]);xlim([0 10]);
% hold on; plot(T,a4,'s','color',[0.1 0.7 0.5])
% ylim([-0.5 1.5])
% figure; plot(T,tmp3*1000,'linewidth',2,'color',[0.1 0.8 0.8])
% yline(b3*1000,'linewidth',1.5)
% hold on; plot(T,b4,'color',[0 0.8 0.6])
% set(gcf,'position',[45 45 1100 200]);xlim([0 10]);ylim([-0.1 1])
% figure; plot(T,c4,'linewidth',2,'color',[0 0.8 0.6]);
% set(gcf,'position',[45 45 1100 200]);xlim([0 10]);ylim([-0.5 1.5])
% set(gca,'xtick',[],'ytick',[])
% [0 0.8 0.6]
% [0.7 0.9 0.5] % yello
%% stage
stg=zeros(length(T),1);
stg(baseline_t*60*fs+1:end)=1;
stg(find(c4==1))=2;
stg(1:baseline_t*60*fs)=0; % erase baseline detection

%figure; plot(T,stg,'k','linewidth',2)

% burst: gap between stage 3 is less than 40 sec
clear tmp* 
tmp=find(bsmtx==1);
tmp1=find(diff(tmp)>fs*40)  % 40sec thresholod
if ~isempty(tmp1)
    for ind=1:length(tmp1)
        if ind==1
            stg(tmp(1):tmp(tmp1(ind)))=3;
        elseif ind~=1
            stg(tmp(tmp1(ind-1)+1):tmp(ind))=3;
        end
        
        if ind==length(tmp1)
            stg(tmp(tmp1(ind)+1):tmp(end))=3;
        end
    end
end
stg(tmp(1):tmp(end))=3;
if (length(stg)-tmp(end))<40*fs %40 second to end of signal length
    stg(tmp(end):end)=3;
end


% gap between stage 2/3 is less than 10 sec -> 2
s1=find(stg>=2);
s2=find(diff(find(stg>=2))>1);
s3=find((s1(s2+1)-s1(s2))<fs*gapsec); % 10sec thresholod
s4=[s1(s2(s3)),s1(s2(s3)+1)]
for ind=1:size(s4,1)
    stg(s4(ind,1):s4(ind,2))=2;
end

%  stage 2 is less than 10 sec -> 1
s1=find(stg==2);
s2=find(diff(find(stg==2))>1);
s3=find((s1([s2;length(s1)] )-s1([1; s2+1]))<fs*gapsec);
s21=[1; s2+1];
s22=[s2;length(s1)];
s4=[s1(s21(s3)),s1(s22(s3))]
for ind=1:size(s4,1)
    stg(s4(ind,1):s4(ind,2))=1;
end


figure; plot(T,stg,'k','linewidth',2)
set(gcf,'position',[45 45 1100 200]);xlim([0 10]);ylim([-0.2 4])
for ind=1:size(sup_datapoint,1)
hold on
plot([T(sup_datapoint(ind,1)) T(sup_datapoint(ind,2))],[3.4 3.4],'r','linewidth',2)
end
end
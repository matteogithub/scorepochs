clc
clear
Frange=1:1:30;
t_win=9760;
fs=160;
n_chan=64;
ep_size=5; %seconds
eps=ep_size*fs;
n_eps=floor((t_win/fs)/ep_size);

addpath /Users/matteofraschini/Downloads/eeglab13_6_5b
inDir='/Users/matteofraschini/Downloads/eegmmidb/';
fil='S*';
cases=dir(fullfile(inDir,fil));

m_score_R01=zeros(length(cases),n_eps);
m_score_R02=zeros(length(cases),n_eps);
B_R01=zeros(length(cases),n_eps);
I_R01=zeros(length(cases),n_eps);
B_R02=zeros(length(cases),n_eps);
I_R02=zeros(length(cases),n_eps);

for i=1:length(cases)
    i
    tic
    EEGR01=pop_biosig(strcat(inDir,cases(i).name,'/',cases(i).name,'R01.edf'), 'importevent','off','importannot','off');
    EEGR02=pop_biosig(strcat(inDir,cases(i).name,'/',cases(i).name,'R02.edf'), 'importevent','off','importannot','off');
    EEGR01=pop_select(EEGR01,'nochannel',{'Status'});
    EEGR02=pop_select(EEGR02,'nochannel',{'Status'});
    
    score_R01=zeros(size(EEGR01.data,1),n_eps);
    score_R02=zeros(size(EEGR02.data,1),n_eps);
    
    for k=1:size(EEGR01.data,1)
        my_data_R01=zeros(n_eps,eps);
        my_data_R02=zeros(n_eps,eps);
        PSD_R01=zeros(n_eps,length(Frange));
        PSD_R02=zeros(n_eps,length(Frange));
        for w=1:n_eps
            end_ep=w*eps;
            in_ep=end_ep-eps+1;
            my_data_R01(w,:)=EEGR01.data(k,in_ep:end_ep);
            my_data_R02(w,:)=EEGR02.data(k,in_ep:end_ep);
            [Pxx_R01,F]=pwelch(my_data_R01(w,:)',[],[],Frange,fs);
            PSD_R01(w,:)=Pxx_R01;
            [Pxx_R02,F]=pwelch(my_data_R02(w,:)',[],[],Frange,fs);
            PSD_R02(w,:)=Pxx_R02;
        end
        c_R01=corr(PSD_R01','type','Spearman');
        c_R01(1:size(c_R01,1)+1:end)=0;
        score_R01(k,:)=sum(c_R01,1)/(size(c_R01,1)-1);
        c_R02=corr(PSD_R02','type','Spearman');
        c_R02(1:size(c_R02,1)+1:end)=0;
        score_R02(k,:)=sum(c_R02,1)/(size(c_R02,1)-1);
    end
    m_score_R01(i,:)=mean(score_R01,1);
    m_score_R02(i,:)=mean(score_R02,1);
    [B_R01(i,:),I_R01(i,:)]=sort(m_score_R01(i,:),'descend');
    [B_R02(i,:),I_R02(i,:)]=sort(m_score_R02(i,:),'descend');
    toc    
end



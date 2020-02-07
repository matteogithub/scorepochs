clc
clear
Frange=1:1:30;
t_win=9760;
fs=160;
n_chan=64;
ep_size=5; %seconds
eps=ep_size*fs;
n_eps=floor((t_win/fs)/ep_size);
lf=1;
hf=30;
alpha=[8 13];

addpath /Users/matteofraschini/Downloads/eeglab13_6_5b
inDir='/Users/matteofraschini/Downloads/eegmmidb/';
fil='S*';
cases=dir(fullfile(inDir,fil));

mean_alphaP_R01=zeros(length(cases),n_eps);
mean_alphaP_R02=zeros(length(cases),n_eps);

for i=1:length(cases)
    i
    tic
    EEGR01=pop_biosig(strcat(inDir,cases(i).name,'/',cases(i).name,'R01.edf'), 'importevent','off','importannot','off');
    EEGR02=pop_biosig(strcat(inDir,cases(i).name,'/',cases(i).name,'R02.edf'), 'importevent','off','importannot','off');
    EEGR01=pop_select(EEGR01,'nochannel',{'Status'});
    EEGR02=pop_select(EEGR02,'nochannel',{'Status'});
    
    for w=1:n_eps
        end_ep=w*eps;
        in_ep=end_ep-eps+1;
        my_data_R01=EEGR01.data(:,in_ep:end_ep);
        my_data_R02=EEGR02.data(:,in_ep:end_ep);
        [Pxx_R01,F]=pwelch(my_data_R01',[],[],Frange,fs);
        PSD_R01=Pxx_R01;        
        [Pxx_R02,F]=pwelch(my_data_R02',[],[],Frange,fs);
        PSD_R02=Pxx_R02;        
        [v,indmin]=min(abs(F-lf));
        [v,indmax]=min(abs(F-hf));        
        [v,indminalpha]=min(abs(F-alpha(1)));
        [v,indmaxalpha]=min(abs(F-alpha(2)));        
        totalP_R01=sum(PSD_R01(indmin:indmax,:));
        alphaP_R01=sum(PSD_R01(indminalpha:indmaxalpha,:))./totalP_R01;
        mean_alphaP_R01(i,w)=mean(alphaP_R01,2);
        totalP_R02=sum(PSD_R02(indmin:indmax,:));
        alphaP_R02=sum(PSD_R02(indminalpha:indmaxalpha,:))./totalP_R02;
        mean_alphaP_R02(i,w)=mean(alphaP_R02,2);
    end
    toc    
end
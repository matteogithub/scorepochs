epl=5;
F=1:1:30;
fs=160;
nep=floor((size(data,2)/fs)/epl);
eps=epl*fs;
score=zeros(size(data,1),nep);
for i=1:size(data,1)
    my_data=zeros(nep,eps);
    PSD=zeros(nep,length(F));
    for j=1:nep
        end_ep=j*eps;
        in_ep=end_ep-eps+1;
        my_data(j,:)=data(i,in_ep:end_ep);
        [Pxx,F]=pwelch(my_data(j,:)',[],[],F,fs);
        PSD(j,:)=Pxx;
    end    
    c=abs(corr(PSD','type','Spearman'));
    c(1:size(c,1)+1:end)=0;
    score(i,:)=sum(c,1)/(size(c,1)-1);   
end
m_score=mean(score,1);
sort(m_score)
plot(m_score,'o');
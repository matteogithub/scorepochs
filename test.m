fs=160;
epl=20;
nep=floor((size(data,2)/fs)/epl);
eps=epl*fs;

score=zeros(size(data,1),nep);
for i=1:size(data,1)
    my_data=zeros(nep,eps);
    for j=1:nep
        end_ep=j*eps;
        in_ep=end_ep-eps+1;
        my_data(j,:)=data(i,in_ep:end_ep);       
    end
    
    c=abs(corr(my_data','type','Spearman'));
    score(i,:)=sum(c,1)/(size(c,1)-1);   
end
m_score=mean(score,1);
sort(m_score)
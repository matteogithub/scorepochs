n_best=4;

I_R01=cell2mat(idx_best_R01);

R01=mean(mean_alphaP_R01(:,I_R01(:,1:n_best)),2);
R02=mean(mean_alphaP_R02(:,I_R02(:,1:n_best)),2);

[P,H,STATS]=signrank(R01,R02);

Cind=combnk([1:1:12],n_best);

dH=zeros(size(Cind,1),1);
dP=zeros(size(Cind,1),1);
dSTATS=zeros(size(Cind,1),1);
for i=1:size(Cind,1)
    i
    cR01=mean(mean_alphaP_R01(:,Cind(i,:)),2);
    cR02=mean(mean_alphaP_R02(:,Cind(i,:)),2);
    [cP,cH,cSTATS]=signrank(cR01,cR02);
    dH(i)=cH;
    dP(i)=cP;
    dSTATS(i)=cSTATS.zval;
end
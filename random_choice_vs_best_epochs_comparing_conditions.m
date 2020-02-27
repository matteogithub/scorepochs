
%     Copyright (C) 2020 Matteo Demuru, Matteo Fraschini
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


fileIN = fullfile(pwd,'output','idxBest_epochs_scores_per_condition.mat');

load(fileIN) 

% compute relative power in alpha band for all the epochs

nSubj = numel(epoch_R01);
nEp   = size(epoch_R01{1},2);

cfg = [];
cfg.freqRange    = 1 : 30; % frequency range of interest to compute the power spectrum
cfg.fs           = 160;    % sample frequency
cfg.windowL      = 5;      % epoch length in second used to segment the data
cfg.smoothFactor = 3;      % smoothing the power spectrum
cfg.freqBOI      = [8 13]; % frequency band of interest for which calculate the relative power

g_rPow_R01       = zeros(nSubj,nEp);
g_rPow_R02       = zeros(nSubj,nEp);

for i = 1 : nSubj 
    
    [~,g_rPow_R01(i,:)] = compute_rel_power(cfg,epoch_R01{i});
    [~,g_rPow_R02(i,:)] = compute_rel_power(cfg,epoch_R02{i});  
end


% compare epochs chosen with scorEpoch with a random selection
% in terms of effect-size

nEp2sel = 4; % subset of epochs to select

idx_cmb = combnk(1:nEp,nEp2sel);

dH     = zeros(size(idx_cmb,1),1);
dP     = zeros(size(idx_cmb,1),1);
dSTATS = zeros(size(idx_cmb,1),1);


figure;
plot(0,0)
hold

max_R02 = zeros(nSubj,1);
min_R01 = ones(nSubj,1);

dist_R01 = zeros(nSubj,size(idx_cmb,1));
dist_R02 = zeros(nSubj,size(idx_cmb,1));

for i = 1 : size(idx_cmb,1)
      
    dist_R01(:,i) = mean(g_rPow_R01(:,idx_cmb(i,:)),2);
    cR01 = mean(g_rPow_R01(:,idx_cmb(i,:)),2);
    
    plot(cR01,'ro');
    
    dist_R02(:,i) = mean(g_rPow_R02(:,idx_cmb(i,:)),2);
    cR02 = mean(g_rPow_R02(:,idx_cmb(i,:)),2);
    
    plot(cR02,'yo');
    
    max_R02 = max(max_R02,cR02);
    min_R01 = min(min_R01,cR01);
    
    
    [cP,cH,cSTATS] = signrank(cR01,cR02);
    
    dH(i)     = cH;
    dP(i)     = cP;
    dSTATS(i) = cSTATS.signedrank;
    
end



idx_chosen_R01 = reshape(cell2mat(idx_best_R01),nEp,nSubj)';
idx_chosen_R02 = reshape(cell2mat(idx_best_R02),nEp,nSubj)';

idx_chosen_R01 = idx_chosen_R01(:,1:nEp2sel);
idx_chosen_R02 = idx_chosen_R02(:,1:nEp2sel);

R01 = mean(g_rPow_R01(:,idx_chosen_R01),2);
R02 = mean(g_rPow_R02(:,idx_chosen_R02),2);

[P,H,STATS] = signrank(R01,R02);

plot(R01,'b*')
plot(R02,'g*')

plot(min_R01,'b+')
plot(max_R02,'g+')

ylabel('Relative alpha','FontSize',16)
xlabel('Subject','FontSize',16)

%% plotting 

% look for normality with Kolmogorov-Smirnov
hR01 = zeros(nSubj,1);
pR01 = zeros(nSubj,1);

hR02 = zeros(nSubj,1);
pR02 = zeros(nSubj,1);

for i = 1 : nSubj 
   
    hR01(i) = kstest(zscore(dist_R01(i,:))); 
    hR02(i) = kstest(zscore(dist_R02(i,:))); 

end

sprintf('R01 Number of subjects with normal distribution %i out of %i ',sum(~hR01),nSubj)
sprintf('R02 Number of subjects with normal distribution %i out of %i ',sum(~hR02),nSubj)

% for each subject check to what percentile the value of the chosen epoch belong to

[sdist_R01,~]     = sort(dist_R01,2,'ascend');
find_chosen_R01   = abs(sdist_R01 - repmat(R01,1,size(sdist_R01,2)));
[~,idx_R01]       = min(find_chosen_R01,[],2);
prcR01            = (idx_R01/size(sdist_R01,2))*100;

[sdist_R02,~]     = sort(dist_R02,2,'ascend');
find_chosen_R02   = abs(sdist_R02 - repmat(R02,1,size(sdist_R02,2)));
[~,idx_R02]       = min(find_chosen_R02,[],2);
prcR02            = (idx_R02/size(sdist_R02,2))*100;



% plot one subject as example

fig(1) = figure;
% selected subject

subj_idx = 1;
nBins   = 20;


[N1 binC1] = hist(dist_R01(subj_idx,:),nBins);
[N2 binC2] = hist(dist_R02(subj_idx,:),nBins);

max_occurences = max([N1 N2]);
max_occurences = max_occurences + 2;

bar(binC1,N1,'FaceColor','r')
hold
bar(binC2,N2,'FaceColor','y')


line([R01(subj_idx) , R01(subj_idx)],[0 , max_occurences],'LineWidth',2);

line([R02(subj_idx) , R02(subj_idx)],[0 , max_occurences],'LineWidth',2);

xlabel('mean relative alpha power','FontSize',14);
ylabel('# Occurences','FontSize',14);


% Plot the overall two distributions across subject (eyes-closed, eyes-open) using all the
% possible combination epochs (i.e. cmb(4,12)). Overlay the two
% distribution using the 'best selected' epochs for each patients


[N1 binC1] = hist(dist_R01(:),nBins);
[N2 binC2] = hist(dist_R02(:),nBins);


[best_N1 best_binC1] = hist(R01,nBins);
[best_N2 best_binC2] = hist(R02,nBins);

max_occurrences_global = max([N1 N2]);
max_occurences_best    = max([best_N1 best_N2]);
x_max                  = 0.7;

fig(2) = figure;

subplot(2,1,1)

bar(binC1,N1/max_occurrences_global,'FaceColor','r','FaceAlpha',0.3)
hold
bar(binC2,N2/max_occurrences_global,'FaceColor','y','FaceAlpha',0.3)
xlim([0 x_max])
title('All possible Combination')
ylabel('Normalized Occurences','FontSize',14);
subplot(2,1,2)

bar(best_binC1,best_N1/max_occurences_best,'FaceColor','b','FaceAlpha',0.3)
hold
bar(best_binC2,best_N2/max_occurences_best,'FaceColor','g','FaceAlpha',0.3)
xlim([0 x_max])
title('Selection using scorEpochs')

xlabel('mean relative alpha power','FontSize',14);
ylabel('Normalized Occurences','FontSize',14);

fig(3) = figure;
plot(binC1,N1/max_occurrences_global,'r')
hold
plot(binC2,N2/max_occurrences_global,'y')

plot(best_binC1,best_N1/max_occurences_best,'b')

plot(best_binC2,best_N2/max_occurences_best,'g')

xlabel('mean relative alpha power','FontSize',14);
ylabel('Normalized Occurences','FontSize',14);
legend({'Eyes-open all combination','Eyes-closed all combination','Eyes-open selection using scorEpochs','Eyes-closed selection using scorEpochs' })











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
idx_cmb = idx_cmb(end:-1:1,:); 

dH          = zeros(size(idx_cmb,1),1);
dP          = zeros(size(idx_cmb,1),1);
dSTATS      = zeros(size(idx_cmb,1),1);
effect_size = zeros(size(idx_cmb,1),1);


dist_R01 = zeros(nSubj,size(idx_cmb,1));
dist_R02 = zeros(nSubj,size(idx_cmb,1));

for i = 1 : size(idx_cmb,1)
      
    dist_R01(:,i) = mean(g_rPow_R01(:,idx_cmb(i,:)),2);
    cR01 = mean(g_rPow_R01(:,idx_cmb(i,:)),2);
    
    
    dist_R02(:,i) = mean(g_rPow_R02(:,idx_cmb(i,:)),2);
    cR02 = mean(g_rPow_R02(:,idx_cmb(i,:)),2);
    

    [cH,cP,~,cSTATS] = ttest(cR01,cR02);
    
    dH(i)          = cH;
    dP(i)          = cP;
    dSTATS(i)      = cSTATS.tstat;
    effect_size(i) = abs(cSTATS.tstat/sqrt(length(cR01)));
    
end

% selection using scorEpochs

idx_chosen_R01 = reshape(cell2mat(idx_best_R01),nEp,nSubj)';
idx_chosen_R02 = reshape(cell2mat(idx_best_R02),nEp,nSubj)';

idx_chosen_R01 = idx_chosen_R01(:,1:nEp2sel);
idx_chosen_R02 = idx_chosen_R02(:,1:nEp2sel);

 R01 = zeros(nSubj,1);
 R02 = zeros(nSubj,1);
for i = 1 : nSubj  
    R01(i,1) = mean(g_rPow_R01(i,idx_chosen_R01(i,:)));
    R02(i,1) = mean(g_rPow_R02(i,idx_chosen_R02(i,:)));
end

[H,P,~,STATS] = ttest(R01,R02);

scorEpochEffSize = abs(STATS.tstat/sqrt(length(R01)));

fig(1) = figure;
subplot(2,3,[1,2])
plot(effect_size,'*--');
hold
line([0 size(idx_cmb,1)],[scorEpochEffSize scorEpochEffSize],'LineStyle','--','Color','green','LineWidth',2.5)
set(gca,'FontSize',30)
xlabel('t-test run','FontSize',33)
ylabel('Cohen''s d effect-size','FontSize',33)
legend({'random selection','scorEpochs selection'},'FontSize',18);
title('(a) Cohen''s d timecourse','FontSize',22)
%fig(2) = figure;
subplot(2,3,[4,5])
[counts,BinC] = hist(effect_size,50);
bar(BinC,counts);
hold
line([scorEpochEffSize scorEpochEffSize],[0 max(counts)+1],'Color','green','LineStyle','--','LineWidth',2.5)
set(gca,'FontSize',30)
xlabel('Cohen''s d','FontSize',33)
ylabel('# Occurrences','FontSize',33)
title('(b) Cohen''s d distribution','FontSize',22);
legend({'random selection','scorEpochs selection'},'FontSize',18)

outFname = {'cohenD_timeline_and_distribution'};
for i = 1 : numel(fig)
  set(fig(i),'WindowState','fullscreen')
  outfile = fullfile(pwd,'output',outFname{i});
  print(fig(i),outfile,'-djpeg');
  close(fig(i))
  
end




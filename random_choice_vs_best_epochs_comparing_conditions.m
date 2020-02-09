
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
for i = 1 : size(idx_cmb,1)
      
    cR01 = mean(g_rPow_R01(:,idx_cmb(i,:)),2);
    
    plot(cR01,'r*');
    cR02 = mean(g_rPow_R02(:,idx_cmb(i,:)),2);
    
    plot(cR02,'y*');
    
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

plot(R01,'bo')
plot(R02,'go')


ylabel('Relative alpha','FontSize',16)
xlabel('Subject','FontSize',16)





% Scorepochs 
%
% Function to select the best (most homogenoous) M/EEG epochs from a
% resting-state recordings. 
%
% INPUT
%    cfg struct with the following fields
%           freqRange    - array with the frequency range used to compute the power
%                          spectrum (see MATLAB pwelch function)
%           fs           - integer representing sample frequency         
%           windowL      - integer representing the window length (in seconds)  
%           smoothFactor - smoothing factor for the power spectrum
% 
%    data             - 2d array with the time-series (channels X time samples)
%     
%
%
% OUTPUT
%      
%    epoch       -  cell array of the data divided in equal length epochs 
%                   of length windowL (channels X time samples)
%                  
%    idx_best_ep - array of indexes sorted according to the best score
%                  this array should be used for the selection of the best
%                  epochs
%
%    score_Xep   - array of score per epoch


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

function [idx_best_ep,epoch,score_Xep] = scorEpochs(cfg,data) 


    % divide the data in epochs of windowL length
    
    epLen   = cfg.windowL * cfg.fs;
    dataLen = size(data,2);
    nCh     = size(data,1);
    idx_ep  = 1 : epLen : dataLen - epLen + 1;
    nEp     = numel(idx_ep);
    
    epoch   = cell(1,nEp);
     
    pxx     = cell(1,nEp);
    
    for e = 1 : nEp 
        
        epoch{e} = data(:,idx_ep(e):idx_ep(e)+epLen-1);
        % compute power spectrum
        pxx{e} = pwelch(epoch{e}',[],[],cfg.freqRange,cfg.fs)';
       
        
        if(cfg.smoothFactor ~= 0)
                pxx{e} = movmean(pxx{e}',cfg.smoothFactor)';
        end
            
    end
    
    % compute score across channels and across epochs 
    pxxXch      = zeros(nEp,numel(cfg.freqRange));
    score_chXep = zeros(nCh,nEp);
    for c = 1 : nCh
        
        for e = 1 : nEp
            
            pxxXch(e,:) = pxx{e}(c,:);  
            
        end
        score_ch = corr(pxxXch','type','Spearman');
        
        score_chXep(c,:) = mean(score_ch);
       
    end    
    
    score_Xep = mean(score_chXep,1);
    
    [~,idx_best_ep] = sort(score_Xep,'descend');
    
    

    
    
    
  
    
    
    

% Compute relative power spectrum (using matlab pwelch for power spectrum)
%
% INPUT
%  
%       
%    cfg struct with the following fields
%           freqRange    - array with the frequency range used to compute the power
%                          spectrum (see MATLAB pwelch function)
%           fs           - integer representing sample frequency         
%           windowL      - integer representing the window length (in seconds)  
%           smoothFactor - smoothing factor for the power spectrum
%           freqBOI      - frequency band of interest for which the
%                          relative power will be calculated
%
%           
% 
%    epoch               - cell array with 2d array consisting of time-series (channels X time samples)
% 
%
%
% OUTPUT
%    
%    
%    rel_pow   - 2d array of relative power spectrum (# channels X # epochs)
%    g_rel_pow - 2d array of global relative across channels (1 X # epochs) 
% 

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


function [rel_pow,g_rel_pow] = compute_rel_power(cfg,epoch)  

boi     = cfg.freqBOI;
nEp     = size(epoch,2);
nCh     = size(epoch{1},1);
rel_pow = zeros(nCh,nEp);

for e = 1 : nEp 
           
    % compute power spectrum
    [pxx , F] = pwelch(epoch{e}',[],[],cfg.freqRange,cfg.fs);
     pxx      = pxx';
    % smoothing the power spectrum
    if(cfg.smoothFactor ~= 0)
        pxx = movmean(pxx',cfg.smoothFactor)';
    end
     
     % compute relative power spectrum
    [~,idx_start] = min(abs(F-boi(1)));
    [~,idx_stop]  = min(abs(F-boi(end)));
   
    rel_pow(:,e) = sum(pxx(:,idx_start:idx_stop),2)./sum(pxx,2);
                 
end

g_rel_pow = mean(rel_pow,1);


% Practical example for scorEpoch using an open-source EEG dataset:
% EEG Motor Movement/Imagery Dataset 
% https://physionet.org/content/eegmmidb/1.0.0/
% Schalk, G., McFarland, D.J., Hinterberger, T., Birbaumer, N., Wolpaw, J.R. BCI2000: A General-Purpose Brain-Computer Interface (BCI) System. 
% IEEE Transactions on Biomedical Engineering 51(6):1034-1043, 2004.

%Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, Mietus JE, Moody GB, Peng C-K, Stanley HE. 
%PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals (2003). Circulation. 101(23):e215-e220.



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


% For each subject resting-state conditions are in R01 and R02 edf files 
% Eyes-open     R01
% Eyes-closed   R02

% Required toolbox EEGLAB to read EDFs
eeglab_path = '/Users/matteofraschini/Downloads/eeglab13_6_5b';

addpath(eeglab_path)
% input folder containing the data from the EEG dataset
inDir       = '/Users/matteofraschini/Downloads/eegmmidb/';   

% scorEpoch input struct

cfg.freqRange    = 1 : 30; % frequency range of interest to compute the power spectrum
cfg.fs           = 160;    % sample frequency
cfg.windowL      = 5;      % epoch length in second used to segment the data
cfg.smoothFactor = 0;      % smoothing factor for the power spectrum

filter   = 'S*';
cases    = dir(fullfile(inDir,filter));
nSubj    = numel(cases); 


idx_best_R01 = cell(1,nSubj); % indexes best epochs ordered (i.e. the best index is in first position)
epoch_R01    = cell(1,nSubj); 
score_R01    = cell(1,nSubj);

idx_best_R02 = cell(1,nSubj);
epoch_R02    = cell(1,nSubj);
score_R02    = cell(1,nSubj);


for i = 1 : length(cases)
    
    
    R01_F  =  fullfile(inDir,cases(i).name,strcat(cases(i).name,'R01.edf')); % file name eyes open
    R02_F  =  fullfile(inDir,cases(i).name,strcat(cases(i).name,'R02.edf')); % file name eyes closed 
    
    % load data (channels X time samples)
    EEGR01 = pop_biosig(R01_F, 'importevent','off','importannot','off');
    EEGR02 = pop_biosig(R02_F, 'importevent','off','importannot','off');
  
    EEGR01 = pop_select(EEGR01,'nochannel',{'Status'});
    EEGR02 = pop_select(EEGR02,'nochannel',{'Status'});
    
    data_R01 = EEGR01.data;
    data_R02 = EEGR02.data;
    
    [idx_best_R01{i},epoch_R01{i},score_R01{i}] = scorEpochs(cfg,data_R01); 
    [idx_best_R02{i},epoch_R02{i},score_R02{i}] = scorEpochs(cfg,data_R02);
  
end



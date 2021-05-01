
function [ score_table ] = pipeline01_dev2_happe(eeglab_dir, code_dir, data_dir, subj_name_edf)
%function [ score_table, eeg_ICA ] = pipeline01_dev2_happe(eeglab_dir, code_dir, data_dir, subj_name_edf)

%    This function performs a series of preprocessing steps
%    comparing two alternative automated pipelines:

%    HAPPE: https://github.com/lcnhappe/happe/blob/master/HAPPE_pipeline_v1_0.m
%    ASR (without data segment removal):
%        https://github.com/sccn/clean_rawdata/wiki

%    and computing the "scopepochs" after each step of the preprocessing

% INPUTS: 
%   subj_name_edf : file with a single subject data in .edf format
%                    (i.e.: subj_name_edf = 'S003R01.edf'
%   optional (already set in the function):
%       eeglab_dir
%            !!! eeglab versione 2021226 requires this list of PLUGIN: 
%                Biosig3.7.5;
%                Cleanline1.04; 
%                clean_rawdata2.3
%                PrepPipeline0.55.4
%                ICLabel
% 
%
%       code_dir (with scorepochs package)
%       subj_dir (also for saving intermediate steps)

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% EXTRA PARAMETERs section to check/set

% OUTPUTs: 
%   score_struct (with score_Xep at each step)

%   andrea.vitale@gmail.com 
%   20210430

%%
    %clear; close all

    % SET DIR
    if ~exist('eeglab_dir') %isempty('eeglab_dir')
        eeglab_dir = 'D:\IIT\EEG\eeglab_20201226'
        cd(eeglab_dir)
        %eeglab
        eeglab('nogui')
    end
    
    if ~exist('code_dir')
        code_dir = 'D:\IIT\_PROJECT\SCORE_epoch\code';
        addpath(genpath(code_dir))
    end

    if ~exist('data_dir')
        data_dir = 'D:\IIT\_PROJECT\SCORE_epoch\data';
        % this folder should contain also the CHANNEL INFO (.txt file)
        cd(data_dir)
    end
    
    if ~exist('subj_name_edf') 
        disp('!!! subj_name is required as INPUT') 
    end
   
    
    %% EXTRA PARAMETERs: - - - - - - - - - - - - - - 
    do_chan_pruning = 1
    do_cleanraw = 1
    do_waveletICA = 1
    
    do_plot_chan = 0
    do_plot_PSD = 1
    
    do_save_cleanline = 0
    do_save_wavclean_ICA = 1
    do_save_cleanraw_avgref_ICA = 1
    do_save_score = 1
    
    
    %% INPUT: resting state EYES OPEN (R01)
    
    % R01 = resting state eyes open - - - - - - - - - - - 
    %subj_name_edf = 'S003R01.edf';
    %file_name = 'S001R01.edf'; 
    %file_name = 'S002R01.edf'; 
    %file_name = 'S003R01.edf'; % score epoch > 95% already for raw data 
    %file_name = 'S010R01.edf'; % score epoch > 90% already for raw data 

    % R02 = resting state eyes close - - - - - - - - - - - 
    %file_name = 'S001R02.edf';

    chan_file = 'coord_BS_motorEEG.txt';

    
%0) LOAD .edf file - - - - - - -
    eeg_struct = []; 
    eeg_struct = pop_biosig(fullfile(data_dir, subj_name_edf));
    %EEG = pop_biosig('D:\IIT\_PROJECT\SCORE_epoch\data\S001R01.edf');
    % eeg_struct = EEG;

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Nose direction should be set from '-Y' to default +X in EEG.chanlocs
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    eeg_struct=pop_chanedit(eeg_struct, 'load',{fullfile(data_dir, chan_file),'filetype','xyz'},'nosedir','-Y');
    %eeg_struct=pop_chanedit(eeg_struct, 'load',{'D:\\IIT\\_PROJECT\\SCORE_epoch\\data\\coord_BS_motorEEG.txt','filetype','xyz'},'nosedir','-Y');
    %EEG=pop_chanedit(EEG, 'lookup','D:\\IIT\\EEG\\eeglab_20201226\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp','load',{'D:\\IIT\\_PROJECT\\SCORE_epoch\\data\\coord_BS_motorEEG.txt','filetype','xyz'},'nosedir','-Y');

    % check PLOT
    if do_plot_chan
        %figure;
        %pop_eegplot( eeg_struct, 1, 1, 1);

        figure; 
        % by NUMBER
        subplot 121
        topoplot([], eeg_struct.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', eeg_struct.chaninfo);
        % by LABEL
        subplot 122
        topoplot([],eeg_struct.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', eeg_struct.chaninfo);
    end
    
    %%
    % SOME CHECKS - - - - - -
    sample_rate = eeg_struct.srate;
    % length in sec of the recording:
    n_sample = eeg_struct.pnts;
    n_sample / sample_rate;  %in sec
    n_chan = eeg_struct.nbchan;

    % number of channel that can be retained for ICA
    %(number of channel)^2 x 20 to 30 data points to perform ICA
    if n_sample > n_chan^2 * 20
        disp([ num2str ' channels can be given as input to ICA'])
    else
        n_chan_max = sqrt(n_sample/20)
        disp([ 'number of channels for ICA should be reduced to ' num2str(n_chan_max)])
    end
    
    
    %CHANNEL PRUNING (to 19 channels):
    chan_toinclude = [ 62, 63, ...
                   48, 52, 56, 53, 49, ...
                   28, 34, 38, 35, 30, ...
                   17, 13, 10, 14, 18, ...
                       3, 4 ]; 
                   % Fpz=64 , Oz=1 not included !!!
                   % Tp9  and Tp10 (mastoids) not included
 
    
% 1) RAW DATA = = = = = = = = = = = = = = = = = = = = = =
    %% remove segment of the data (at the end of the recording) with 0 values
    % !!!! specific for BCI2000 dataset
    i_chan = 1;
    for i_sample = 1:n_sample
        if eeg_struct.data(i_chan,end-i_sample) ~= 0
            zero_last_sample = n_sample-i_sample;
            break
        end
    end
    
    % reject last timepoints
    eeg_struct = eeg_eegrej( eeg_struct, [zero_last_sample n_sample] );
    eeg_raw = eeg_struct;
    
    
    % = = = = = = = = = = = = = = = = = = = = = = == = =
    %% RUN SCOREPOCHS at each step of the preprocessing:
    %INPUT
    %    cfg struct with the following fields
    %           freqRange    - array with the frequency range used to compute the power
    %                          spectrum (see MATLAB pwelch function)
    %           fs           - integer representing sample frequency         
    %           windowL      - integer representing the window length (in seconds)  
    %           smoothFactor - smoothing factor for the power spectrum
    cfg = []; 
    % <<<<<<<<<<<<<<<<<< ENTIRE FREQUENCY RANGE<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %cfg.freqRange = [ 1 : 80 ];
    % <<<<<<<<<<<<<<<<<< ONLY ALPHA BAND <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cfg.freqRange = [ 8 : 13 ]; 
    cfg.fs = eeg_struct.srate;
    cfg.windowL = 5; % in sec <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cfg.smoothFactor = 0;

       
% ----------------------------------------------------------------
%2) BAND-PASS FILTERED data 
    
    hpf_cutoff = 1;
    %lpf_cutoff = 80;
    lpf_cutoff = eeg_struct.srate/2 -1;
    
    line_noise_freq = 60; %<<<<<<<<<<< TO SET <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
    eeg_struct = pop_eegfiltnew(eeg_struct, hpf_cutoff, [], [],0,[],0);
    eeg_hpf = eeg_struct;
    
    % - - - 
    eeg_struct = pop_eegfiltnew(eeg_struct, [], lpf_cutoff, [],0,[],0);
    eeg_lpf = eeg_struct;
        
% ----------------------------------------------------------------
    
%     if do_save_cleanline
%         cd(fullfile(data_dir))
%         pop_saveset(eeg_cleanline, 'filename', [subj_name_edf(1:end-4) '_cleanline'])
%     end
    % load
    % eeg_cleanline = EEG;
    
    
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% ----------------------------------------------------------------
% Two alternative sub-pipelines
% ----------------------------------------------------------------
if do_waveletICA
    
%3  CLEANLINE (line noise removal)
    % https://github.com/sccn/cleanline
    % If cleaning continuous (un-epoched data) 
    % then you may wish to use sliding windows of 3-4 seconds with 50% overlap

    eeg_cleanline = pop_cleanline(eeg_struct, 'Bandwidth',2,'ChanCompIndices',[1:eeg_struct.nbchan] , ...
        'SignalType','Channels','ComputeSpectralPower',true,'LineFrequencies',[line_noise_freq] , ...
        'NormalizeSpectrum',false,'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,...
        'ScanForLines',true,'SmoothingFactor',100,'VerbosityLevel',1,'SlidingWinLength',...
        eeg_struct.pnts/eeg_struct.srate,'SlidingWinStep',eeg_struct.pnts/eeg_struct.srate);
    %eeg_cleanline = eeg_struct;
    
%   in HAPPE
%   EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',chan_index,'computepower',1,'linefreqs',...
%         [60 120] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
%         'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');


%4a  WAVELET-ICA
    % run wavelet-ICA (ICA first for clustering the data, then wavelet thresholding on the ICs)
    %uses a soft, global threshold for the wavelets, wavelet family is coiflet (level 5), threshold multiplier .75 to remove more high frequency noise
    %for details, see wICA.m function
    
    EEG = eeg_cleanline;
    
    %% crude bad channel detection using spectrum criteria and 3SDeviations as channel outlier threshold, done twice
    EEG = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],'threshold',[-3 3],'norm','on','measure','spec','freqrange',[hpf_cutoff lpf_cutoff]); %1 125
    %EEG.setname='rawEEG_f_cs_ln_badc';
    
    EEG = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 125]);
    
    % INTERPOLATE:
    EEG = pop_interp(EEG, eeg_cleanline.chanlocs, 'spherical');
       
    %rank(eeg_struct.data)
    EEG = pop_select(EEG, 'channel', chan_toinclude);
    eeg_psdthresh_badchan_interp = EEG;
    
    try 
        %if pipeline_visualizations_semiautomated == 0
        %    [wIC, A, W, IC] = wICA(EEG,'runica', 1, 0, [], 5);
        %elseif pipeline_visualizations_semiautomated == 1
            [wIC, A, W, IC] = wICA(EEG,'runica', 1, 1, EEG.srate, 5);
        %end
    catch wica_err
        if strcmp ('Output argument "wIC" (and maybe others) not assigned during call to "wICA".',wica_err.message)
            error('Error during wICA, most likely due to memory settings. Please confirm your EEGLAB memory settings are set according to the description in the HAPPE ReadMe')
        else
            rethrow(wica_err)
        end
    end
        
    %reconstruct artifact signal as channelsxsamples format from the wavelet coefficients
    artifacts = A*wIC;
    
    %reshape EEG signal from EEGlab format to channelsxsamples format
    EEG2D=reshape(EEG.data, size(EEG.data,1), []);
    
    %subtract out wavelet artifact signal from EEG signal
    EEG.data = EEG2D-artifacts;
    %eeg_wavclean = EEG2D-artifacts;

    eeg_wavclean = EEG; 
    
% ----------------------------------------------------------------
% 5a    
    eeg_wavclean_ICA = pop_runica(eeg_wavclean, 'icatype', 'runica', 'extended',1,'interrupt','on');

    %%(PLOT component topography:)
    % (see: https://github.com/sccn/viewprops)
    %pop_topoplot(eeg_ICA, 0, [1:length(chan_toinclude)] ,'EDF file',[5 5] ,0,'electrodes','on');
    
    eeg_wavclean_ICA = pop_iclabel(eeg_wavclean_ICA, 'default');

    %pop_viewprops(eeg_wavclean_ICA, 0, [1:eeg_wavclean_ICA.nbchan], [2 80]) % for component properties
    
    if do_save_wavclean_ICA
        cd(fullfile(data_dir))
        pop_saveset(eeg_wavclean_ICA, 'filename', [subj_name_edf(1:end-4) '_wavclean_ICA'])
    end
                        
    eeg_wavclean_nobadICA = pop_icflag(eeg_wavclean_ICA, [0 0.2;0.7 1;0.7 1;NaN NaN;0.7 1;0.7 1;0.7 1]);
   

%6a RE-REFERENCE
    eeg_wavclean_nobadICA_avgref = pop_reref(eeg_wavclean_nobadICA, []);
    
    
% !!! TO DO
%     %store IC variables and calculate variance of data that will be kept after IC rejection:
%     ICs_to_keep =find(EEG.reject.gcompreject == 0);
%     ICA_act = EEG.icaact;
%     ICA_winv =EEG.icawinv;
    
    
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
end

%elseif do_cleanraw
if do_cleanraw
    
%3  NOTCH FILTER:
    eeg_notch = pop_eegfiltnew(eeg_struct, 'locutoff',line_noise_freq-2, ...
                              'hicutoff',line_noise_freq+2,'revfilt',1,'plotfreqz',1);


%4b CLEANRAW data
    % methodA (suggested in "HAPPE")
    % crude bad channel detection using PSD +/- 3SDeviations as channel outlier threshold, done twice
    % EEG = pop_rejchan(EEG, 'elec',chan_index,'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 125]);
    
    % or 
    % methodB
    % with CLEANRAW data (and Artifact Subspace Reconstruction ASR)
    % https://github.com/sccn/clean_rawdata/wiki
    
    % on the channel pruned dataset:
    %eeg_struct = eeg_cleanline;
    eeg_cleanraw = pop_clean_rawdata(eeg_notch, ...
                    'FlatlineCriterion',5, ...
                    'ChannelCriterion',0.8, ...
                    'LineNoiseCriterion',4, ...
                    'Highpass','off','WindowCriterion','off', ...
                    'BurstCriterion',20, 'BurstRejection','off','Distance','Euclidian');
                    % !!! burst criterion OFF -> in order to repair but not to remove (cut) datapoint
    %eeg_cleanraw = eeg_struct;
        
    % why this step is better implemented in this way ?? 
    % (see: Comments to the HAPPE paper, and how to choose the critical parameters
    
    %ASR == good at removing occasional large-amplitude noise/artifacts
    %ICA == good at decomposing constant fixed-source noise/artifacts/signals

    % vis_artifacts to compare the cleaned data to the original.

   
% ----------------------------------------------------------------
%5b BAD CHANNEL (DETECTED by cleanraw and) INTERPOLATED
    % !!! NOT CLEAR IF BAD CHANNELs have been already removed (by cleanraw function
       
    % NO SAMPLE rejected 
%     bad_sample_cleanraw = eeg_cleanraw.etc.clean_sample_mask;
%     % percentage of data kept
%     bad_data_percent = sum(bad_sample_cleanraw)/eeg_struct.pnts*100
%     disp(['percentage of data suggested for removal = % ' num2str(bad_data_percent)])
    
    bad_chan_cleanraw = eeg_cleanraw.etc.clean_channel_mask;
    bad_chan_label = {};
    ii = 1;
    for i_chan = 1:length(bad_chan_cleanraw)
        if bad_chan_cleanraw(i_chan) == 0
            %bad_chan_label = [ bad_chan_label, eeg_struct.chanlocs(i_chan).labels ];
            bad_chan_label{ii} = eeg_cleanraw.chanlocs(i_chan).labels;
            bad_chan_idx(ii) = i_chan;
            ii = ii+1;
        end
    end
    
    %cleanraw_badchan = {bad_chan_label}
    
    %if do_interp_badchan
        % INTERPOLATE instead of remove
        % based on the channel location of the initial eeg_struct
        eeg_cleanraw_badchan_interp = pop_interp(eeg_cleanraw, eeg_struct.chanlocs, 'spherical');
        %eeg_cleanraw = pop_interp(eeg_cleanraw, bad_chan_idx, 'spherical');
    %end
    % ?? do NOT INTERPOLATE before ICA ??
    
    
% ----------------------------------------------------------------
%6  RE-REFERENCE to the AVERAGE    
    eeg_cleanraw_avgref = pop_reref(eeg_cleanraw_badchan_interp, []);
    %eeg_avgref = eeg_struct;
    
% ----------------------------------------------------------------
%7  CHAN PRUNING + IC decomposition (and rejection)    
    
    eeg_cleanraw_avgref_red = pop_select(eeg_cleanraw_avgref, 'channel', chan_toinclude);
    

    eeg_cleanraw_avgref_ICA = pop_runica(eeg_cleanraw_avgref_red, 'icatype', 'runica', 'extended',1,'interrupt','on');

    %%(PLOT component topography:)
    % (see: https://github.com/sccn/viewprops)
    %pop_topoplot(eeg_ICA, 0, [1:length(chan_toinclude)] ,'EDF file',[5 5] ,0,'electrodes','on');
    
    eeg_cleanraw_avgref_ICA = pop_iclabel(eeg_cleanraw_avgref_ICA, 'default');
                     % for component viewing
    %pop_viewprops(eeg_cleanraw_ICA, 0, [1:eeg_cleanraw_ICA.nbchan], [2 50]) % for component properties
    %pop_viewprops( eeg_ICA, 0, [1:length(chan_toinclude)], [2 80], [], 0, eeg_ICA.etc.ic_classification, 'on')
        %spec_opt, erp_opt, scroll_event, classifier_name, fig)

    if do_save_cleanraw_avgref_ICA
        cd(fullfile(data_dir))
        pop_saveset(eeg_cleanraw_avgref_ICA, 'filename', [subj_name_edf(1:end-4) '_cleanraw_avgref_ICA'])
    end
    
    % - - - -  - - - - - - - - - - - - - - - -
    % REMOVE BAD COMPONENT based on ICLabels
    %eeg_nobadICA = pop_icflag(eeg_ICA, [NaN NaN;0.8 1;0.8 1;NaN NaN;0.8 1;0.8 1;0.8 1]);

    % if the % of brain ICA < 0.2 -> then is removed
    eeg_cleanraw_avgref_nobadICA = pop_icflag(eeg_cleanraw_avgref_ICA, [0 0.2;0.7 1;0.7 1;NaN NaN;0.7 1;0.7 1;0.7 1]);
    %eeg_nobadICA = EEG;

    % or manually: 
    % bad_ICA_idx = []; 
    % eeg_struct = pop_subcomp(eeg_struct, bad_ICA_idx, 0);

end

    if do_plot_PSD
        % CHECK the PSD (before and after LINE NOISE removal)                                
        figure; subplot 131; %!!!
        pop_spectopo(eeg_struct, 1, [ ], 'EEG' , 'percent', 50, 'freq', [8 13 20], 'freqrange',[2 lpf_cutoff],'electrodes','off');

        subplot 132;
        pop_spectopo(eeg_cleanline, 1, [ ], 'EEG' , 'percent', 50, 'freq', [8 13 20], 'freqrange',[2 lpf_cutoff],'electrodes','off');
        
        subplot 133;
        pop_spectopo(eeg_notch, 1, [ ], 'EEG' , 'percent', 50, 'freq', [8 13 20], 'freqrange',[2 lpf_cutoff],'electrodes','off');
    end
    
    
    % = ==========================================================
    %% COMPUTE SCOREPOCHS at each step:
    % = ==========================================================
    
    prep_step = {
            'eeg_raw';
            'eeg_hpf';
            'eeg_lpf';
            
            % alternative a) - - - - - - - - -
    
            'eeg_cleanline';
            'eeg_psdthresh_badchan_interp';
            'eeg_wavclean';
            'eeg_wavclean_nobadICA';
            'eeg_wavclean_nobadICA_avgref';
            
            % alternative b) - - - - - - - - -
    
            'eeg_notch';
            'eeg_cleanraw';
            'eeg_cleanraw_badchan_interp';
            'eeg_cleanraw_avgref';
            'eeg_cleanraw_avgref_nobadICA';
                 };
             
    % CREATE A SCORE STRUCT for final report:
    % with epoch not sorted !!!
    score_struct = [];
    
    for i_step = 1:length(prep_step)
        eval(['eeg_step_tmp = ' prep_step{i_step} ]);
        
        % reduce to 19 che number of channels 
        if eeg_step_tmp.nbchan > length(chan_toinclude)
            eeg_step_tmp =  pop_select(eeg_step_tmp, 'channel', chan_toinclude);
        end
        
        [idx_best_ep,epoch,score_Xep] = scorEpochs(cfg, eeg_step_tmp.data);
        eval([ 'score_struct.' prep_step{i_step} '= score_Xep' ]);
        disp(mean(score_Xep))
    end
  
    
    if do_save_score
        save_name = [subj_name_edf(1:end-4) '_scorepoch' ];
        save(save_name, 'score_struct')
        
%         % REPORT other METRICS:
%         n_chan; 
%         chan_interpolated = []; 
%         n_badICA =[];
%         Percent_Variance_Kept_of_Post_Waveleted_Data=[];

    end
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end
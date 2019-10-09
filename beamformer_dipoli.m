% Fieldtrip - Beamformer
% Script written by: Sara Jamil
% 
% EEG data info:
% Control subject - mTBI_P19 
% P300 ERPs
% Standard tone - S10 

display('Would you like to clear all and run the whole beamformer_dipoli program?');
pause;

clear all; close all; clc;

%% Import data and preprocess

% Import data
fn = 'C:\Users\Sara\Documents\Fieldtrip\beamforming_test_(sample_data)\mTBI_P19.bdf';
cfg = [];
cfg.dataset = fn;
cfg.trialfun = 'ft_trialfun_general';


cfg.trialdef.eventtype = 'STATUS';
cfg.trialdef.prestim = 0.1;
cfg.trialdef.poststim = 0.9;
cfg.trialdef.eventvalue = [10, 1034]; 
%10,1034-control; 20,1044-deviant?; 30,1054-duration deviant;
%40,1064-deviant?
cfg = ft_definetrial(cfg);

% Select trials to remove MANUALLY
trials = 493:1:2462; % P300: 10,1034 - 492 trials
% trials = 421:1:1940; % P300: 10 - 420 trials
cfg.trl(trials,:) = []; 

% Preprocessing
cfg.reref = 'yes';
cfg.refchannel = {'all'}; %'EXG2'. 'EXG3' or 'all' for common avg ref *****
cfg.hpfilter = 'yes'; % high-pass filter 0.5Hz
cfg.hpfreq = 0.5;
cfg.lpfilter = 'yes'; % low-pass filter 30Hz
cfg.lpfreq = 30;
cfg.lpfiltord = 10; % increase to get a more ideal filter 
cfg.demean = 'yes';
data = ft_preprocessing(cfg);

% Add matching 10-20 naming scheme
elecs = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1'...
    ,'C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7'...
    ,'PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','AFz','Fz'...
    ,'F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8'...
    ,'TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
data.label(1:64) = elecs';

% Use only eeg channels
cfg = [];
cfg.channel = 'eeg';
data = ft_preprocessing(cfg, data);

% % Save the preprocessed, trial split data
save data data

% % % Load data
% load data


%% Covariance calculation

%%% Timelockanalysis
% Calculating the covariance matrix
% For common spatial filter
cfg = [];
cfg.covariance = 'yes';
cfg.channel = 'eeg';
cfg.vartrllength = 0; %default=0 (trials all same length)
cfg.covariancewindow = 'all';
cfg.keeptrials = 'yes';
tlock = ft_timelockanalysis(cfg, data);

% % Save the timelock
save tlock tlock

% % % Load the tlock
% load tlock


%% Preparing the MRI

% % % Segment the template MRI
% % mri = ft_read_mri('C:\Users\Sara\Documents\Fieldtrip\fieldtrip-20160104\template\anatomy\single_subj_T1.nii');
% % mri.coordsys = 'spm';
mri = ft_read_mri('Subject01.mri');
% cfg = [];
% cfg.output = {'scalp', 'skull', 'brain'};
% mri_seg = ft_volumesegment(cfg, mri);
% 
% % % Save the segmented data
% save mri_seg mri_seg

% % Load saved volume segmentation
load mri_seg


% % % Mesh the MRI for bemcp or dipoli Headmodel
% cfg = [];
% cfg.tissue = {'scalp','skull','brain'};
% cfg.method = 'projectmesh';
% cfg.numvertices = [1000 2000 3000];
% mri_mesh = ft_prepare_mesh(cfg, mri_seg);
% 
% % % Save the mesh
% save mri_mesh mri_mesh

% % Load the mesh
% load mri_mesh


%% Preparing the Headmodel

% Head model
% 'bemcp'
% cfg = [];
% cfg.method = 'bemcp';
% cfg.tissue = {'scalp','skull','brain'};
% cfg.conductivity = [1 1/80 1]*0.33; %scalp:skull:brain conductivities
% vol = ft_prepare_headmodel(cfg, mri_mesh);


% % Cannot use 'dipoli' method on windows

% Load the dipoli volume conduction model
load mri_vol_dipoli.mat
% 
% Convert units to 'mm'
vol = ft_convert_units(vol, 'mm');


%% Aligning the electrodes

% % Electrodes
% elec = ft_read_sens('standard_1020.elc');
% elec = ft_convert_units(elec, 'mm'); % Same units as MRI
% 
% % Plot the geometry of the volume conduction model
% figure;
% ft_plot_mesh(vol.bnd(1), 'edgecolor', 'none', 'facealpha', 0.8, 'facecolor', [0.6 0.6 0.8]);
% hold on;
% ft_plot_mesh(vol.bnd(2), 'edgecolor', 'none', 'facealpha', 0.4);
% hold on;
% ft_plot_mesh(vol.bnd(3), 'facecolor', [0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05); 
% hold on;
% ft_plot_sens(elec, 'style', 'sk');
% hold on;
% camlight;
% 
% % Align the electrodes
% nas = mri.hdr.fiducial.mri.nas;
% lpa = mri.hdr.fiducial.mri.lpa;
% rpa = mri.hdr.fiducial.mri.rpa;
% 
% trans = mri.transform;
% nas = ft_warp_apply(trans, nas, 'homogenous');
% lpa = ft_warp_apply(trans, lpa, 'homogenous');
% rpa = ft_warp_apply(trans, rpa, 'homogenous');
% 
% % Create a structure similar to a template set of electrodes
% fid.chanpos = [nas; lpa; rpa];
% fid.label = {'Nz', 'LPA', 'RPA'};
% fid.unit = 'mm';
% 
% % Alignment
% cfg = [];
% cfg.method = 'fiducial';
% cfg.template = fid;
% cfg.elec = elec;
% cfg.fiducial = {'Nz', 'LPA', 'RPA'};
% elec_align = ft_electroderealign(cfg);
% 
% % % Save the electrode alignment
% % save elec_align elec_align
% 
% % % Load the electrode alignment
% % load elec_align
% 
% % Check the alignment again
% figure;
% ft_plot_mesh(vol.bnd(1), 'edgecolor', 'none', 'facealpha', 0.8, 'facecolor', [0.6 0.6 0.8]);
% hold on;
% ft_plot_mesh(vol.bnd(2), 'edgecolor', 'none', 'facealpha', 0.4);
% hold on;
% ft_plot_mesh(vol.bnd(3), 'facecolor', [0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05); 
% hold on;
% ft_plot_sens(elec, 'style', 'sk');
% hold on;
% camlight;
% 
% 
% % Realign electores - interactive
% % Rotate [23 0 4]
% % Scale [1 1 1]
% % Translate [15 17 0]
% cfg = [];
% cfg.method = 'interactive';
% cfg.elec = elec_align;
% cfg.headshape = vol.bnd(1);
% elec_realign = ft_electroderealign(cfg);
% % NEW:
% % Rotate [0 0 0]
% % Scale [1 1 1]
% % Translate [15 0 -5]
% 
% % % Save the realigned electrodes
% save elec_realign elec_realign

% % Load realigned electrodes
load elec_realign

% % Plot the REALIGNED ELECTRODES
% % 'dipoli'
% figure;
% ft_plot_mesh(vol.bnd(1), 'edgecolor', 'none', 'facealpha', 0.8, 'facecolor', [0.6 0.6 0.8]);
% hold on;
% ft_plot_mesh(vol.bnd(2), 'edgecolor', 'none', 'facealpha', 0.4);
% hold on;
% ft_plot_mesh(vol.bnd(3), 'facecolor', [0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05); 
% hold on;
% ft_plot_sens(elec_realign, 'style', 'sk');
% hold on;
% camlight;



%% Preparing the leadfield matrix

%%% Compute the Leadfield
cfg = [];
cfg.elec = elec_realign;
cfg.senstype = 'eeg';
cfg.headmodel = vol;
cfg.channel = {'eeg'};
cfg.grid.resolution = 5; % use a 3D grid with a 5 mm resolution
cfg.grid.unit = 'mm';
cfg.grid.tight = 'yes';
cfg.reducerank = 3; % default = 3 for EEG, removes the weakest orientation
cfg.normalize = 'no'; % default = 'no'; otherwise need normalizeparam (0.5)
% if you are not contrasting the activity of interest against another
% condition or baseline time-window
cfg.backproject = 'yes'; % default = 'yes'; 
% determines when reducerank is applied whether the lower rank leadfield
% is projected back to the original linear subspace or not
[grid] = ft_prepare_leadfield(cfg, tlock); % tlock contains trials
% Sourcemodel is computed with leadfield

% % Save the leadfield
save grid grid

% % Load leadfield
% load grid



%% Beamforming

%%% LCMV Beamformer
% Compute common spatial filter
cfg = [];
cfg.elec = elec_realign;
cfg.senstype = 'eeg';
cfg.channel = {'eeg'};
cfg.method = 'lcmv';
cfg.headmodel = vol; % headmodel
cfg.grid = grid; % leadfield
cfg.keeptrials = 'yes';
cfg.keepleadfield = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.projectnoise = 'yes'; % noise estimate computed on the bases of the 
% % smallest eigenvalue of the CSD (or covariance)
% % Used to remove center of head bias
cfg.rawtrial = 'no'; % Construct filter from single trials, apply to
% % single trials. Note: keeptrials if using grid.filter
source_lcmv = ft_sourceanalysis(cfg, tlock);

% % Save the common filter
save source_lcmv source_lcmv

% % Load the common filter
% load source_lcmv




% % % % Try using sourcedescriptives
% % % % Result: fail
% % % % cannot trace if source_lcmv powmethod is not trace
% % % % powmethod trace is disabled due to bugzilla bug 2395 ?
% % cfg = [];
% % cfg.powmethod = 'trace';
% % cfg.projectmom = 'no';
% % source_desc = ft_sourcedescriptives(cfg, source_lcmv);

% % % % Try using sourcemovie
% % % % Result: fail
% % % % Error: source.tri missing, 
% % % % this function requires a triangulated cortical sheet as source model
% % cfg = [];
% % cfg.funparameter = 'pow';
% % cfg.maskparameter = cfg.funparameter;
% % ft_sourcemovie(cfg, source_lcmv);

% % % Project trials through common spatial filter
% % % WILL CALCULATE AVERAGE POWER OVER ALL TIME
% % % TAKES TOO LONG TO CALCULATE ALL TRIALS
% cfg = [];
% cfg.elec = elec_realign;
% cfg.senstype = 'eeg';
% cfg.channel = {'eeg'};
% cfg.method = 'lcmv';
% cfg.headmodel = vol;
% cfg.grid = grid;
% cfg.grid.filter = source_lcmv.avg.filter;
% % % use common spatial filter computed over entire trial length
% cfg.keeptrials = 'yes';
% cfg.rawtrial = 'yes'; 
% % % construct filter from single trials, apply to single trials. 
% % % Note that you also may want to set cfg.keeptrials='yes' to keep all trial information, 
% % % especially if using in combination with grid.filter
% cfg.lcmv.keepfilter = 'yes';
% cfg.lcmv.projectnoise = 'yes';
% source_lcmv_trials = ft_sourceanalysis(cfg, tlock);
% % 
% % % % Save the trial sources
% % save source_lcmv_trials source_lcmv_trials
% % 
% % % % Load
% % % load source_lcmv_trials
% % 

%% Time Windows
% Split trials into time windows
% Can take average over all trials for each time window
% Then perform source analysis with source_lcmv filter

%%% Redefine trials
% Split each trial into 16 windows
% of length 32 samples
% without any overlap
for t = 1:length(data.trialinfo)
    t
    cfg = [];
    cfg.trials = t;
    cfg.length = 0.0625; %32/512
    cfg.overlap = 0;
    datawin16{t} = ft_redefinetrial(cfg, data);
end

%% Timelockanalysis on windows
% Calculating the covariance matrix
% for all windows
% Average of each time window is taken
for w = 1:length(datawin16{1,t}.trialinfo) % 16
    w
    datatemp = data;
    % replace trial and time
    for i = 1:length(datatemp.trialinfo) % for all 492 trials
        datatemp.trial{1,i} = datawin16{1,i}.trial{1,w}; % for window 1 {1,w}
        datatemp.time{1,i} = datawin16{1,i}.time{1,w};
        % all trials collected by window
    end
    cfg = 0;
    cfg.trials = 'all';
    cfg.covariance = 'yes';
    cfg.channel = 'eeg';
    cfg.vartrllength = 0; %default=0 (trials all same length)
    cfg.covariancewindow = 'all';
    %cfg.keeptrials = 'yes'; %no need to keep trials 
    tlockwin16{w} = ft_timelockanalysis(cfg, datatemp);
end

% % Save the timelock data
save tlockwin16 tlockwin16

% % % Load the timelock data
% load tlockwin16


%% Leadfield calculation on windows
% % % Unnecessary to calculate leadfield for each window
% % % Should be the same
% % 
% % for w = 1:length(tlockwin16)
% %     w
% %     cfg = [];
% %     cfg.elec = elec_realign;
% %     cfg.senstype = 'eeg';
% %     cfg.headmodel = vol;
% %     cfg.channel = {'eeg'};
% %     cfg.grid.resolution = 5; % use a 3D grid with a 5 mm resolution
% %     cfg.grid.unit = 'mm';
% %     cfg.grid.tight = 'yes';
% %     cfg.reducerank = 3; % default = 3 for eeg removes the weakest orientation
% %     cfg.normalize = 'no'; % default = 'no'; otherwise need normalizeparam (0.5)
% %     cfg.backproject = 'yes'; % default = 'yes'; 
% %     gridwin16{w} = ft_prepare_leadfield(cfg, tlockwin16{1,w});
% % end
% % 
% % % % % Save the leadfield
% % % save gridwin16 gridwin16
% % 
% % % % Load the leadfield
% % load gridwin16


%% Beamform windows
% % Source analysis using source_lcmv average filter
% % apply to averaged time windows
% % UPDATE: Use leadfield for whole trial length
for w = 1:length(tlockwin16)
    w
    cfg = [];
    cfg.elec = elec_realign;
    cfg.senstype = 'eeg';
    cfg.channel = {'eeg'};
    cfg.method = 'lcmv';
    cfg.headmodel = vol; % headmodel
    cfg.grid = grid; % leadfield
    cfg.grid.filter = source_lcmv.avg.filter;
    %cfg.keeptrials = 'yes';
    %cfg.keepleadfield = 'yes';
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes'; % noise estimate computer on the bases of the 
    % % smallest eigenvalue of the CSD (or covariance?)
    % % Used to remove center of head bias
    cfg.rawtrial = 'no'; % Construct filter from single trials, apply to
    % % single trials. Note: keeptrials if using grid.filter
    sourcelcmvwin16{w} = ft_sourceanalysis(cfg, tlockwin16{1,w});
end

% % Save the sourcewin
display('Time it took to save sourcelcmvwin16:');
tic;
save sourcelcmvwin16 sourcelcmvwin16
toc;

% % % Load the sourcewin
% load sourcelcmvwin16


%% Neural activity index for Windows
% % anatomy originally uses mri_slice instead of mri (faster?)
% % % mri_slice = ft_volumereslice([], mri);
for w = 1:length(sourcelcmvwin16)
    w
    % computing sources without center of head bias
    sourceNAIwin16{1,w} = sourcelcmvwin16{1,w};
    sourceNAIwin16{1,w}.avg.pow = ...
        sourcelcmvwin16{1,w}.avg.pow ./ sourcelcmvwin16{1,w}.avg.noise;
    % plot the source
    cfg = [];
    cfg.parameter = 'pow';
    interpwin16{1,w} = ft_sourceinterpolate(cfg, sourceNAIwin16{1,w}, mri);
end

% save all the windows separately because the variable is too big
% each file has the variable saved as 'temp'
for i = 1:length(interpwin16)
    display(['Time it took to save interpwin',int2str(i),':']);
    tic;
    filename = ['interpwin16_',int2str(i),'.mat'];
    temp = interpwin16{1,i};
    save(filename, 'temp');
    toc;
end

% % load the interp data
% for i = 1:16 % must input this manually...
%     filename = ['interpwin16_',int2str(i),'.mat'];
%     load(filename);
%     interpwin16{1,w} = temp;
%     clear temp
% end


%% Plotting the slices
for i = 1:length(interpwin16)
    i
    name = ['plotslicewin16_',int2str(i)];
    cfg = [];
    cfg.method = 'slice';
    % cfg.anaparameter = mri;
    cfg.funparameter = 'pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.funcolorlim = [0.0 65.0]; % 'zeromax' or range [5.0 25.0]
    cfg.opacitylim = [0.0 65.0];
    cfg.opacitymap = 'rampup';
    ft_sourceplot(cfg, interpwin16{1,i});
    print(name, '-dpng');
end

%%% Putting it into a movie
for i = 1:length(interpwin16)
    name = ['plotslicewin16_',int2str(i),'.png'];
    img = imread(name);
    frame(i) = im2frame(img);
end
% save as .avi file
% choose the frame rate
movie2avi(frame, 'moviewin16_2fps', 'fps', 2);
movie2avi(frame, 'moviewin16_4fps', 'fps', 4);
movie2avi(frame, 'moviewin16_6fps', 'fps', 6);
movie2avi(frame, 'moviewin16_8fps', 'fps', 8);
movie2avi(frame, 'moviewin16_10fps', 'fps', 10);
movie2avi(frame, 'moviewin16_12fps', 'fps', 12);
movie2avi(frame, 'moviewin16_14fps', 'fps', 14);
movie2avi(frame, 'moviewin16_15fps', 'fps', 15);
movie2avi(frame, 'moviewin16_16fps', 'fps', 16);


% Plot reconstruct spatial map for both the main task and the control task
% Author: Zoe Duan 
% Date: 04/21/23

clear; clc;

% add paths
addpath(genpath('/System/Volumes/Data/d/DATA/home/zoe/vistasoft_ts/'));
addpath(genpath('/System/Volumes/Data/d/DATA/home/zoe/preproc_mFiles/'));
addpath('/Users/zoe/Github/princeton-mvpa-toolbox'); 
mvpa_add_paths;

% define expname
expname = 'WM';
% expname = 'control';

% define subjects
subjIDs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
n = length(subjIDs);

% define paths
mainpath = input('Input the path for the parent folder', 's');
% mainpath = '/System/Volumes/Data/d/DATC/datc/formatwm_ZD/public/';

if strcmp(expname, 'WM') 
    datapath_fmri = [mainpath, 'data_fmri_WMTask/'];
    datapath_MVPA = [mainpath, 'results_decoding_WMTask/'];
    newpath = [mainpath, 'results_spaceRecon_WMTask/'];
    % define epochs
    epoch = 'DELAY';
elseif strcmp(expname, 'control') 
    datapath_fmri = [mainpath, 'data_fmri_controlTask/'];
    datapath_MVPA = [mainpath, 'results_decoding_controlTask/'];
    newpath = [mainpath, 'results_spaceRecon_controlTask/'];
    % define epochs
    epoch = 'STIM';
end


% define hemisphere you want to analysis
HEMIs = {'bilat'};     % HEMIs = ['bilat', 'lh', 'rh'];
% define ROIs you want to draw
ROIs = {'V1', 'V2', 'V3', 'V3AB', 'IPS0_IPS1', 'IPS2_IPS3', 'sPCS', 'iPCS'};

% define parameters for spatial reconstruction
% voxel selection parameters
para.ecclim = 20;  % limits for eccentricity
para.sigmalim = 100; % sigmalim=100: no selection
para.velim = 0; % velim=0: no selection based on variance explained
% reconstruction map parameters
para.binunit = 0.1; % downsampling
para.sigmafix = 5;
% plotting recon parameters
para.stimr = 6; % radius of stim circle, consistant with task design
% threshold for fitting line
para.perc = 10;
% line width for stimulus and fitting line
para.lineWidth = 5.0;

% define modulators
moduls = {'angular', 'radial'};


%% combine all subjects data together if n>1
if n > 1
    % initializa variables
    recon_angular_allSubjs = [];
    recon_radial_allSubjs = [];

    % Cycle over subj
    for sii = 1:n
        subjID = subjIDs(sii);

        % define path for loading data
        datapath = [newpath, 'subj', num2str(subjID), '/'];

        % load recon for weights from classifier
        load([datapath, 'recon_angular_', epoch, '.mat']);
        load([datapath, 'recon_radial_', epoch, '.mat']);

        if isempty(recon_angular_allSubjs)
            recon_angular_allSubjs = recon_angular;
            recon_radial_allSubjs = recon_radial;
        else
            recon_angular_allSubjs = cat(3, recon_angular_allSubjs, recon_angular);
            recon_radial_allSubjs = cat(3, recon_radial_allSubjs, recon_radial);
        end
    end % of cyclying over subj

    % save combined results for all subjects
    outpath = [newpath, 'allSubjs/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % save recon results into disk
    save([outpath, 'recon_angular_', epoch, '.mat'], 'recon_angular_allSubjs', '-v7.3');
    save([outpath, 'recon_radial_', epoch, '.mat'], 'recon_radial_allSubjs', '-v7.3');

end



%% ploting spatial reconstruction maps

% define figure formats
% fig_form = 'separate';       % separate each orientation condition
fig_form = 'combine';       % rotate all orientation conditions to align them centered at zero

% define line-fitting form
% fit_form = '';
% fit_form = 'fitline_all';   % fit lines based on all pixels on the image
fit_form = 'fitline_stim';  % fit lines only based on pixels within the stimulus

% define path to save figures
if strcmp(fig_form, 'separate')
    para.mapSize = [0 0 1200 1600];
    outpath_fig = [newpath, 'Figures_separate/'];
elseif strcmp(fig_form, 'combine')
    para.mapSize = [0 0 700 1600];
    outpath_fig = [newpath, 'Figures_combine/'];
end
if ~exist(outpath_fig, 'dir')
   mkdir(outpath_fig)
end

% define datapath for loading data
if n == 1
    datapath = [newpath, 'subj', num2str(subjIDs(1)), '/'];
else
    datapath = [newpath, 'allSubjs/'];
end
% load recon for real data
load([datapath, 'recon_angular_', epoch, '.mat']);
load([datapath, 'recon_radial_', epoch, '.mat']);

% Cycle over hemispheres
for hii = 1:length(HEMIs)
    hemi = HEMIs{hii};

    % Cycle over ROIs
    for rii = 1:length(ROIs)
        roi_name = ROIs{rii};

        % create recon maps for current ROI
        recon = cell(1, length(moduls), n);

        % plot spatial map for current ROI
        if n == 1
            % delete empty [] in recon,
            % which doesn't have data for this ROI for this subj
            emptyCells = all(cellfun(@isempty, recon_angular), 1);
            recon_angular = recon_angular(:, ~emptyCells);
            emptyCells = all(cellfun(@isempty, recon_radial), 1);
            recon_radial = recon_radial(:, ~emptyCells);

            % create recon maps for current ROI
            recon(1,1,1) = recon_radial(hii,rii);
            recon(1,2,1) = recon_angular(hii,rii);

            % save subject-wised plots
            final_fig = [outpath_fig, 'subj', num2str(subjIDs(1)), '/'];
            if ~exist(final_fig, 'dir')
                mkdir(final_fig)
            end
            savePath = [final_fig, hemi, '_', roi_name, '_spatialMaps_', epoch, '_real_subj', num2str(subjIDs(1))];
        else
            recon(1,1,:) = recon_radial_allSubjs(hii,rii,:);
            recon(1,2,:) = recon_angular_allSubjs(hii,rii,:);

            % save allSubjs plots
            savePath = [outpath_fig, hemi, '_', roi_name, '_spatialMaps_', epoch, '_real'];
        end

        if strcmp(fig_form, 'separate')
            plot_recon_sep(recon, para, [], rii, roi_name, savePath, fit_form);
        elseif strcmp(fig_form, 'combine')
            plot_recon_com(recon, para, [], rii, roi_name, savePath, fit_form);
        end

    end % of cyclying ROIs
end % of cyclying hemis
close all



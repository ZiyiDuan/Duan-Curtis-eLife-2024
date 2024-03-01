% Orientation decoding analysis for cross tasks
% Data source is percentage signal change(PSC) or beta coefficents
% Analysis is based on Princeton MVPA toolbox (http://www.csbmb.princeton.edu/mvpa).
% Please download the toolbox and add them to the paths before running this code.
% Analysis including:
%   - Train classifiers from control task and test classifiers on WM task
%   - Within vs. across modulators
% Author: Zoe Duan
% Date: 05/15/23

clear; clc;

% add paths
restoredefaultpath
addpath('/Users/zoe/Github/princeton-mvpa-toolbox');
mvpa_add_paths;

% define parameters
subjIDs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
n = length(subjIDs);

% define hemisphere you want to analysis
HEMIs = {'bilat'};     % HEMIs = ['bilat', 'lh', 'rh'];
% define ROIs you want to draw
ROIs = {'V1', 'V2', 'V3', 'V3AB', 'IPS0_IPS1', 'IPS2_IPS3', 'sPCS', 'iPCS'};

% define paths
mainpath = input('Input the path for the parent folder', 's');
% mainpath = '/System/Volumes/Data/d/DATC/datc/formatwm_ZD/public/';

datapath_fmri = [mainpath, 'data_fmri_WMTask/'];
datapath_fmri_control = [mainpath, 'data_fmri_controlTask/'];
newpath = [mainpath, 'results_decoding_crossTask/'];
if ~exist(newpath, 'dir')
    mkdir(newpath)
end

% define classifier
class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.act_funct{1} = 'softmax'; % 'logsig';
class_args.performFcn = 'crossentropy'; % 'msesparse';
class_args.nHidden = 0;

% define # of iterations for cross dataset classification
nIterations = 10;


%% MVPA decoding cross tasks
% Cycle over subj
for sii = 1:n
    subjID = subjIDs(sii);
    sprintf(['Decoding for subj', num2str(subjID), '...'])

    % define output path for current subj
    outpath = [newpath, 'subj', num2str(subjID), '/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};

            % define rootpath
            rootpath = [datapath_fmri, 'subj', num2str(subjID)];
            rootpath_control = [datapath_fmri_control, 'subj', num2str(subjID)];
            
            % decide if data exist
            filename = [rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'];
            if ~exist(filename, 'file')
                sprintf([hemi, '_', roi_name, ' does not exist'])
            else
                sprintf(['Decoding for ', hemi, '_', roi_name])

                %% load data structure into disk
                % load data for main task
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'])
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_radial.mat'])
                subj_angular_WM = subj_angular;
                subj_radial_WM = subj_radial;
                % load data for control task
                load([rootpath_control, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'])
                load([rootpath_control, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_radial.mat'])
                subj_angular_control = subj_angular;
                subj_radial_control = subj_radial;

                %% load voxel inds for different task
                % load inds for WM task
                load([rootpath, '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_voxel_inds_STIM.mat'])
                load([rootpath, '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_voxel_inds_DELAY.mat'])
                voxel_inds_WM_DELAY = voxel_inds_DELAY;
                % load inds for control task
                load([rootpath_control, '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_voxel_inds_STIM.mat'])
                voxel_inds_control_STIM = voxel_inds_STIM;
                
                % find the shared inds for WM task and control task
                [~, shared_inds_WM, shared_inds_control] = intersect(voxel_inds_WM_DELAY, voxel_inds_control_STIM, 'rows');

                %% decoding within modulator
                % for angular modulator
                [results_DELAY_angular] = decoding_crossTask(subj_angular_control, 'STIM', 'oriRef', 'runs', shared_inds_control, ...
                    subj_angular_WM, 'DELAY', 'oriRef', 'runs', shared_inds_WM, class_args, nIterations);
                save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_DELAY_angular.mat'], 'results_DELAY_angular', '-v7.3')
                % for radial modulator
                [results_DELAY_radial] = decoding_crossTask(subj_radial_control, 'STIM', 'oriRef', 'runs', shared_inds_control, ...
                    subj_radial_WM, 'DELAY', 'oriRef', 'runs', shared_inds_WM, class_args, nIterations);
                save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_DELAY_radial.mat'], 'results_DELAY_radial', '-v7.3')

                %% decoding across modulator
                % for a2r modulator
                [results_DELAY_a2r] = decoding_crossTask(subj_angular_control, 'STIM', 'oriRef', 'runs', shared_inds_control, ...
                    subj_radial_WM, 'DELAY', 'oriRef', 'runs', shared_inds_WM, class_args, nIterations);
                save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_DELAY_a2r.mat'], 'results_DELAY_a2r', '-v7.3')
                % for r2a modulator
                [results_DELAY_r2a] = decoding_crossTask(subj_radial_control, 'STIM', 'oriRef', 'runs', shared_inds_control, ...
                    subj_angular_WM, 'DELAY', 'oriRef', 'runs', shared_inds_WM, class_args, nIterations);
                save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_DELAY_r2a.mat'], 'results_DELAY_r2a', '-v7.3')

            end
        end % of cycling over ROIs
    end % of cyclying over hemis
end % of cyclying subjects





% Reconstruct spatial map for both tasks
% refer to code here: https://github.com/clayspacelab/abstract_format_WM
% Author: Zoe Duan 
% Date: 10/16/23


clear; clc;

% add paths
addpath('/Users/zoe/Github/princeton-mvpa-toolbox'); 
mvpa_add_paths;

% define expname
expname = 'WM';
% expname = 'control';

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

datapath_pRFs = [mainpath, 'data_pRFs/'];
if strcmp(expname, 'WM') 
    datapath_fmri = [mainpath, 'data_fmri_WMTask/'];
    newpath = [mainpath, 'results_spaceRecon_WMTask/'];
    % define epochs
    epoch = 'DELAY';
elseif strcmp(expname, 'control') 
    datapath_fmri = [mainpath, 'data_fmri_controlTask/'];
    newpath = [mainpath, 'results_spaceRecon_controlTask/'];
    % define epochs
    epoch = 'STIM'; 
end
if ~exist(newpath, 'dir')
   mkdir(newpath)
end

% define parameters for spatial reconstruction
% voxel selection parameters
para.ecclim = 20;  % limits for eccentricity
para.sigmalim = 100; % sigmalim=100: no selection
para.velim = 0; % velim=0: no selection based on variance explained
% reconstruction map parameters
para.binunit = 0.1; % downsampling
para.sigmafix = 5; % 2.5, 5, 7.5
% plotting recon parameters
para.stimr = 6; % radius of stim circle, consistant with task design

%% Spatial reconstruction 
% Cycle over subj
for sii = 1:n
    subjID = subjIDs(sii);
    sprintf(['Reconstructing for subj', num2str(subjID), '...'])

    % initialize recon variables for each subjects and ROIs
    recon_angular = cell(length(HEMIs), length(ROIs));
    recon_radial = cell(length(HEMIs), length(ROIs));

    % load pRF parameters
    pRFpath = [datapath_pRFs, 'subj', num2str(subjID), '/pRFfit/'];
    % get pRF parameters for current subj
    pRFparams = niftiread([pRFpath, 'RF_ss5_25mm-fFit.nii.gz']);

    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};
            roi_filename = [hemi, '.', roi_name, '.nii.gz'];

            % define rootpath
            rootpath_MVPA = [datapath_fmri, 'subj', num2str(subjID)];
            
            % decide if data exist
            filename = [rootpath_MVPA, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'];
            if ~exist(filename, 'file')
                sprintf([hemi, '_', roi_name, ' does not exist'])
            else
                % load voxel indices for current ROI
                load([rootpath_MVPA, '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_voxel_inds_', epoch, '.mat'])
                if strcmp(epoch, 'STIM')
                    voxel_inds = voxel_inds_STIM;
                elseif strcmp(epoch, 'DELAY')
                    voxel_inds = voxel_inds_DELAY;
                end
                
                % load data structure into disk
                load([rootpath_MVPA, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'])
                load([rootpath_MVPA, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_radial.mat'])
                
                % do reconstruction
                % for angular modulator
                recon = do_recon_real(subj_angular, epoch,'oriRef','runs',pRFparams,voxel_inds,'lim',para);
                recon_angular{hii, rii} = recon;
                % for radial modulator
                recon = do_recon_real(subj_radial, epoch,'oriRef','runs',pRFparams,voxel_inds,'lim',para);
                recon_radial{hii, rii} = recon;

            end
        end % of cyclying ROIs
    end % of cyclying hemis

    % save recon results into disk
    % define output path for current subj
    outpath = [newpath, 'subj', num2str(subjID), '/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end
    save([outpath, 'recon_angular_', epoch, '.mat'], 'recon_angular', '-v7.3');
    save([outpath, 'recon_radial_', epoch, '.mat'], 'recon_radial', '-v7.3');

end % of cyclying subjects


% Reconstruct spatial map for model simulations of the two modulated 
% gratings, lines, and dots
% refer to code here: https://github.com/clayspacelab/abstract_format_WM
% Author: Zoe Duan
% Date: 08/13/23

clear; clc;

% add paths
addpath('/Users/zoe/Github/princeton-mvpa-toolbox'); 
mvpa_add_paths;

% define the system
exptSystem = 'scanRoom';

% define model output source
model_source = 'maxLev';
% model_source = 'avg';

% define paths
mainpath = input('Input the path for the parent folder', 's');
% mainpath = '/System/Volumes/Data/d/DATC/datc/formatwm_ZD/public/';
datapath = [mainpath, 'modelResults/'];
if strcmp(model_source, 'maxLev')
    newpath = [datapath, '/results_spaceRecon_model_maxLev/'];
elseif strcmp(model_source, 'avg')
    newpath = [datapath, '/results_spaceRecon_model_avg/'];
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

% load maxLev information
load([datapath, '/modulDiff_maxLev.mat']);




%% Spatial reconstruction for angular modulator
% load data saved by simulate_voxelResp.m
load([datapath, '/voxelResp_angular.mat']);

% Cycle over subj
for sii = 1:numel(subjIDs)
    subjID = subjIDs(sii);
    sprintf(['Reconstructing for subj', num2str(subjID), '...'])

    % initialize recon variables for each subjects and ROIs
    recon_angular = cell(length(HEMIs), length(ROIs));

    % define pRFpath for current subj
    pRFpath = [mainpath, 'pRFs/subj', num2str(subjID), '/pRFfit/'];
    % get pRF parameters for current subj
    pRFparams = niftiread([pRFpath, 'RF_ss5_25mm-fFit.nii.gz']);

    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};
            roi_filename = [hemi, '.', roi_name, '.nii.gz'];

            % find the model response for current subject
            % Only choose oris 15, 75, 135
            if strcmp(model_source, 'maxLev') 
                model_resp = squeeze(rfRespLev_angular{sii,hii,rii}(3:5, maxLev, :));
            elseif strcmp(model_source, 'avg')
                model_resp = squeeze(rfResp_angular{sii,hii,rii}(3:5, :));
            end
            % do spatial reconstruction
            recon = do_recon_model(model_resp, pRFparams,voxel_inds_angular{sii},'lim',para);
            recon_angular{hii,rii} = recon{1};

        end % of cyclying ROIs
    end % of cyclying hemis

    % save recon results into disk
    % define output path for current subj
    outpath = [newpath, 'subj', num2str(subjID), '/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % save recon results into disk
    save([outpath, 'recon_angular.mat'], 'recon_angular', '-v7.3');

end % of cyclying subjects




%% Spatial reconstruction for radial modulator
% load data saved by simulate_voxelResp.m
load([datapath, '/voxelResp_radial.mat']);

% Cycle over subj
for sii = 1:numel(subjIDs)
    subjID = subjIDs(sii);
    sprintf(['Reconstructing for subj', num2str(subjID), '...'])

    % initialize recon variables for each subjects and ROIs
    recon_radial = cell(length(HEMIs), length(ROIs));

    % define pRFpath for current subj
    pRFpath = [mainpath, 'pRFs/subj', num2str(subjID), '/pRFfit/'];
    % get pRF parameters for current subj
    pRFparams = niftiread([pRFpath, 'RF_ss5_25mm-fFit.nii.gz']);

    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};
            roi_filename = [hemi, '.', roi_name, '.nii.gz'];

            % find the model response for current subject
            % Only choose oris 15, 75, 135
            if strcmp(model_source, 'maxLev')
                model_resp = squeeze(rfRespLev_radial{sii,hii,rii}(3:5, maxLev, :));
            elseif strcmp(model_source, 'avg')
                model_resp = squeeze(rfResp_radial{sii,hii,rii}(3:5, :));
            end
            % do spatial reconstruction
            recon = do_recon_model(model_resp, pRFparams,voxel_inds_radial{sii},'lim',para);
            recon_radial{hii,rii} = recon{1};

        end % of cyclying ROIs
    end % of cyclying hemis

    % save recon results into disk
    % define output path for current subj
    outpath = [newpath, 'subj', num2str(subjID), '/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % save recon results into disk
    save([outpath, 'recon_radial.mat'], 'recon_radial', '-v7.3');

end % of cyclying subjects






%% Spatial reconstruction for lines 
% load data saved by simulate_voxelResp.m
load([datapath, '/voxelResp_lines.mat']);

% Cycle over subj
for sii = 1:numel(subjIDs)
    subjID = subjIDs(sii);
    sprintf(['Reconstructing for subj', num2str(subjID), '...'])

    % initialize recon variables for each subjects and ROIs
    recon_lines = cell(length(HEMIs), length(ROIs));

    % define pRFpath for current subj
    pRFpath = [mainpath, 'pRFs/subj', num2str(subjID), '/pRFfit/'];
    % get pRF parameters for current subj
    pRFparams = niftiread([pRFpath, 'RF_ss5_25mm-fFit.nii.gz']);

    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};
            roi_filename = [hemi, '.', roi_name, '.nii.gz'];

            % find the model response for current subject
            % Only choose oris 15, 75, 135
            if strcmp(model_source, 'maxLev') 
                model_resp = squeeze(rfRespLev_lines{sii,hii,rii}(3:5, maxLev, :));
            elseif strcmp(model_source, 'avg')
                model_resp = squeeze(rfResp_lines{sii,hii,rii}(3:5, :));
            end
            % do spatial reconstruction
            recon = do_recon_model(model_resp, pRFparams,voxel_inds_lines{sii},'lim',para);
            recon_lines{hii,rii} = recon{1};

        end % of cyclying ROIs
    end % of cyclying hemis

    % save recon results into disk
    % define output path for current subj
    outpath = [newpath, 'subj', num2str(subjID), '/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % save recon results into disk
    save([outpath, 'recon_lines.mat'], 'recon_lines', '-v7.3');

end % of cyclying subjects





%% Spatial reconstruction for one dot 
% load data saved by simulate_voxelResp.m
load([datapath, '/voxelResp_1dot.mat']);

% Cycle over subj
for sii = 1:numel(subjIDs)
    subjID = subjIDs(sii);
    sprintf(['Reconstructing for subj', num2str(subjID), '...'])

    % initialize recon variables for each subjects and ROIs
    recon_dots = cell(length(HEMIs), length(ROIs));

    % define pRFpath for current subj
    pRFpath = [mainpath, 'pRFs/subj', num2str(subjID), '/pRFfit/'];
    % get pRF parameters for current subj
    pRFparams = niftiread([pRFpath, 'RF_ss5_25mm-fFit.nii.gz']);

    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};
            roi_filename = [hemi, '.', roi_name, '.nii.gz'];

            % find the model response for current subject
            % Only choose oris 15, 75, 135
            if strcmp(model_source, 'maxLev') 
                model_resp = squeeze(rfRespLev_dots{sii,hii,rii}(:, maxLev, :));
            elseif strcmp(model_source, 'avg')
                model_resp = squeeze(rfResp_dots{sii,hii,rii}(:, :));
            end
            % do spatial reconstruction
            recon = do_recon_model(model_resp, pRFparams,voxel_inds_dots{sii},'lim',para);
            recon_dots{hii,rii} = recon{1};

        end % of cyclying ROIs
    end % of cyclying hemis

    % save recon results into disk
    % define output path for current subj
    outpath = [newpath, 'subj', num2str(subjID), '/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % save recon results into disk
    save([outpath, 'recon_1dot.mat'], 'recon_dots', '-v7.3');

end % of cyclying subjects


% Use each voxel's pRF as a mask to sample model output for dots
% for model-based spatial reconstruction
% for both one dot and two dots images
% Author: Zoe Duan
% Date: 02/21/2024


clear; clc

% define the system
exptSystem = 'scanRoom';

% load the screenP parameters based on the system
load(['screenP', '_', exptSystem, '.mat']); 
pix2deg = round(screenP.pixels_per_deg_width);

% define subjs
subjIDs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
n = length(subjIDs);

% define hemisphere you want to analysis
HEMIs = {'bilat'};     % HEMIs = ['bilat', 'lh', 'rh'];
% define ROIs you want to draw
ROIs = {'V1'};

% define paths
mainpath = input('Input the path for the parent folder', 's');
% mainpath = '/System/Volumes/Data/d/DATC/datc/formatwm_ZD/public/';
datapath = [mainpath, 'modelResults/'];

% define parameters for selecting voxels
selection_cutoff = 0.1;     % select voxels with VE>selection_cutoff


%% simulate voxel responses for one dot
% load data saved by generate_modelOutput_dots.m
load([datapath, '/modelOutput_1dot.mat']);

% Create grids for defining a Gaussian pRF function
[Y, X] = meshgrid(1:dims(1),1:dims(2)); % Be careful here is [Y,X] instead of [X,Y]
centerX = (max(X(:))-min(X(:)))/2;
centerY = (max(Y(:))-min(Y(:)))/2;
%convert X Y pixel values to degrees
X = (X-centerX)./pix2deg;
Y = (Y-centerY)./pix2deg;

% change model outputs to vectors
for i = 1:numel(oris)
    vecBands_dots(:,i) = reshape(sumBands_dots(:,:,i),1,dims(1)*dims(2));
    for lev=1:numLevels
        vecBandsLev_dots(:, i, lev) = reshape(sumBandsLev_dots{lev}(:,:,i),1,dims(1)*dims(2));
    end
end

% simulate voxel response based on each subj's pRF parameters
% Cycle over subj
for sii = 1:n
    subjID = subjIDs(sii);
    sprintf(['Simulating for subj', num2str(subjID), '...'])

    % define ROIpath for current subj
    ROIpath = [mainpath, '/pRFs/subj', num2str(subjID), '/ROIs/'];
    % define pRFpath for current subj
    pRFpath = [mainpath, '/pRFs/subj', num2str(subjID), '/pRFfit/'];
    % get pRF parameters for current subj
    pRFparams = niftiread([pRFpath, 'RF_ss5_25mm-fFit.nii.gz']);
    
    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};
            roi_filename = [hemi, '.', roi_name, '.nii.gz'];
            
            % decide whether the ROI exist
            if ~exist([ROIpath, roi_filename], 'file')
                sprintf([ROIpath, roi_filename, ' does not exist'])
            else
                % get the ROI you want to choose
                ROI = niftiRead([ROIpath roi_filename]);
                % get variance explained from pRF fitting
                VE = pRFparams(:, :, :, 2);
                % record selected voxel indices
                voxel_inds_orig = reshape(1:size(VE,1)*size(VE,2)*size(VE,3)*size(VE,4),size(VE));
                % select voxels based on ROI.data
                VE = VE(ROI.data>0);
                voxel_inds_orig = voxel_inds_orig(ROI.data>0);
                % select voxels based on selection cutoff
                voxel_inds = voxel_inds_orig(VE>=selection_cutoff);
                n_voxels = numel(voxel_inds);
                % save voxel inds for each subj and each roi
                voxel_inds_dots{sii,hii,rii} = voxel_inds;

                % get pRF parameters for current ROI based on voxel_inds
                sigma = pRFparams(:, :, :, 4);
                sigma = sigma(voxel_inds);
                x0 = pRFparams(:, :, :, 6);
                x0 = x0(voxel_inds);
                y0 = pRFparams(:, :, :, 7);
                y0 = y0(voxel_inds);

                % calculate model's tuning curve for each voxel pRF
                rfMaskVec = zeros(n_voxels,dims(1)*dims(2),'single');
                for vox=1:n_voxels
                    % define 2D Gaussian pRF mask
                    rfMask = single(exp(-((X-x0(vox)).^2 + (Y-y0(vox)).^2)/(2*sigma(vox)^2))'); 
                    %%%%% why use transpose here???? 
                    % Because the size is (1920, 1080) for rfMask, but the
                    % dims of image is (1080, 1920)
                    rfMaskVec(vox,:) = reshape(rfMask,1,dims(1)*dims(2));
                end

                % initialize voxel responses, (n_oris, n_voxels)
                rfResp_dots{sii,hii,rii} = zeros(numel(oris),n_voxels);
                rfRespLev_dots{sii,hii,rii} = zeros(numel(oris),numLevels,n_voxels);
                % calculate voxel responses, (n_conds, n_voxels)
                rfResp_dots{sii,hii,rii}(:,:) = (rfMaskVec*vecBands_dots)';
                for lev=1:numLevels
                    rfRespLev_dots{sii,hii,rii}(:,lev,:) = (rfMaskVec*vecBandsLev_dots(:, :, lev))';
                end
            end
        end % of cycling over ROIs
    end % of cycling over hemi
end % of cycling over subj

% save results
save([datapath, '/voxelResp_1dot.mat'], 'subjIDs', 'HEMIs', 'ROIs', 'voxel_inds_dots', 'rfResp_dots', 'rfRespLev_dots', 'oris', 'numLevels', '-v7.3');


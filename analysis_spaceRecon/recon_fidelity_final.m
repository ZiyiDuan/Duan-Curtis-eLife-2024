% Compute filtered responses and reconstruction fidelity for both tasks
% Author: Zoe Duan
% Date: 10/17/23

clear; clc;

% add paths
addpath(genpath('/System/Volumes/Data/d/DATA/home/zoe/vistasoft_ts/'));
addpath(genpath('/System/Volumes/Data/d/DATA/home/zoe/preproc_mFiles/'));
addpath('/Users/zoe/Github/princeton-mvpa-toolbox');
mvpa_add_paths;


% define expname
% expname = 'WM';
expname = 'control';

% define line filter form
% filter_form = 'full';     % use all pixels of the recon image
filter_form = 'stim';     % use pixels within the stimulu

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

if strcmp(expname, 'WM') 
    datapath_main = [mainpath, 'results_spaceRecon_WMTask/'];
    % define epochs
    epoch = 'DELAY';
elseif strcmp(expname, 'control')
    datapath_main = [mainpath, 'results_spaceRecon_controlTask/'];
    % define epochs
    epoch = 'STIM'; 
end

% define parameters
ori_labels = [15 75 135]; 
para.ori_filter_width = 1e3;     % projected distance squared limit, used to constrain filter width
para.binunit = 0.1; % downsampling
para.stimr = 6; % radius of stim circle, consistant with task design


%% compute filtered responses and reconstruction fidelity
% Cycle over subj
for sii = 1:n
    subjID = subjIDs(sii);
    sprintf(['Compute recon fidelity for subj', num2str(subjID), '...'])

    % initialize variables
    filterResp_angular = cell(length(HEMIs), length(ROIs));
    filterResp_radial = cell(length(HEMIs), length(ROIs));

    fidelity_angular = cell(length(HEMIs), length(ROIs));
    fidelity_radial = cell(length(HEMIs), length(ROIs));

    % define datapath for loading data
    datapath = [datapath_main, 'subj', num2str(subjID), '/'];

    % load recon for real data
    if strcmp(epoch, 'STIM')
        load([datapath, 'recon_angular_STIM.mat']);
        load([datapath, 'recon_radial_STIM.mat']);
    elseif strcmp(epoch, 'DELAY')
        load([datapath, 'recon_angular_DELAY.mat']);
        load([datapath, 'recon_radial_DELAY.mat']);
    end

    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};

            % find current recon
            recon_angular_curr = recon_angular{hii, rii};
            recon_radial_curr = recon_radial{hii, rii};

            % compute recon fidelity
            [filterResp_angular{hii, rii}, fidelity_angular{hii, rii}] = ...
                compute_fidelity_real(recon_angular_curr, ori_labels, para, filter_form);
            [filterResp_radial{hii, rii}, fidelity_radial{hii, rii}] = ...
                compute_fidelity_real(recon_radial_curr, ori_labels, para, filter_form);

        end % of cyclying ROIs
    end % of cyclying hemis

    % save results into disk
    save([datapath, 'filterResp_angular_', epoch, '_filterForm_', filter_form, '.mat'], 'filterResp_angular', '-v7.3');
    save([datapath, 'filterResp_radial_', epoch, '_filterForm_', filter_form, '.mat'], 'filterResp_radial', '-v7.3');
    save([datapath, 'fidelity_angular_', epoch, '_filterForm_', filter_form, '.mat'], 'fidelity_angular', '-v7.3');
    save([datapath, 'fidelity_radial_', epoch, '_filterForm_', filter_form, '.mat'], 'fidelity_radial', '-v7.3');

end % of cyclying subjects


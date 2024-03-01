% Orientation decoding analysis for the WM and control task separately
% Data source is percentage signal change(PSC) or beta coefficents
% Analysis is based on Princeton MVPA toolbox (http://www.csbmb.princeton.edu/mvpa). 
% Please download the toolbox and add them to the paths before running this code.
% Analysis including:
%   - Decoding orientations for DELAY periods within each modulator
%   - Decoding orientations for DELAY periods cross modulators
%   - Getting weights from classifier by using all data
% Author: Zoe Duan 
% Date: 02/12/23

clear; clc;

% add paths
restoredefaultpath
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

if strcmp(expname, 'WM')
    datapath_fmri = [mainpath, 'data_fmri_WMTask/'];
    newpath = [mainpath, 'results_decoding_WMTask/'];
    % define epochs
    epoch = 'DELAY';
elseif strcmp(expname, 'control')
    datapath_fmri = [mainpath, 'data_fmri_controlTask/'];
    newpath = [mainpath, 'results_decoding_controlTask/'];
    % define epochs
    epoch = 'STIM';
end
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

%% MVPA decoding within modulator

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

            % decide if data exist
            filename = [rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'];
            if ~exist(filename, 'file')
                sprintf([hemi, '_', roi_name, ' does not exist'])
            else
                sprintf(['Decoding for ', hemi, '_', roi_name])
                % load data structure into disk
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'])
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_radial.mat'])

                %% decoding 
                if strcmp(epoch, 'STIM')
                    % for angular modulator
                    [results_STIM_angular] = decoding_within(subj_angular,'STIM','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_STIM_angular.mat'], 'results_STIM_angular', '-v7.3')
                    % for radial modulator
                    [results_STIM_radial] = decoding_within(subj_radial,'STIM','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_STIM_radial.mat'], 'results_STIM_radial', '-v7.3')

                elseif strcmp(epoch, 'DELAY')
                    % for angular modulator
                    [results_DELAY_angular] = decoding_within(subj_angular,'DELAY','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_DELAY_angular.mat'], 'results_DELAY_angular', '-v7.3')
                    % for radial modulator
                    [results_DELAY_radial] = decoding_within(subj_radial,'DELAY','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_DELAY_radial.mat'], 'results_DELAY_radial', '-v7.3')
                end

            end
        end % of cycling over ROIs
    end % of cyclying over hemis
end % of cyclying subjects

                            

%% MVPA decoding cross modulators
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

            % decide if data exist
            filename = [rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'];
            if ~exist(filename, 'file')
                sprintf([hemi, '_', roi_name, ' does not exist'])
            else
                % load data structure into disk
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'])
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_radial.mat'])

                %% decoding 
                if strcmp(epoch, 'STIM')
                    % from angular to radial modulator
                    [results_STIM_a2r] = decoding_cross(subj_angular, 'STIM', 'oriRef', 'runs', ...
                        subj_radial, 'STIM', 'oriRef', 'runs', class_args, nIterations);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_STIM_a2r.mat'], 'results_STIM_a2r', '-v7.3')
                    % from radial to angular modulator
                    [results_STIM_r2a] = decoding_cross(subj_radial, 'STIM', 'oriRef', 'runs', ...
                        subj_angular, 'STIM', 'oriRef', 'runs', class_args, nIterations);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_STIM_r2a.mat'], 'results_STIM_r2a', '-v7.3')
                elseif strcmp(epoch, 'DELAY')
                    % from angular to radial modulator
                    [results_DELAY_a2r] = decoding_cross(subj_angular, 'DELAY', 'oriRef', 'runs', ...
                        subj_radial, 'DELAY', 'oriRef', 'runs', class_args, nIterations);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_DELAY_a2r.mat'], 'results_DELAY_a2r', '-v7.3')
                    % from radial to angular modulator
                    [results_DELAY_r2a] = decoding_cross(subj_radial, 'DELAY', 'oriRef', 'runs', ...
                        subj_angular, 'DELAY', 'oriRef', 'runs', class_args, nIterations);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_DELAY_r2a.mat'], 'results_DELAY_r2a', '-v7.3')
                end

            end
        end % of cycling over ROIs
    end % of cyclying over hemis
end % of cyclying subjects



%% get weights from classification by using all data
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

            % decide if data exist
            filename = [rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'];
            if ~exist(filename, 'file')
                sprintf([hemi, '_', roi_name, ' does not exist'])
            else
                % load data structure into disk
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_angular.mat'])
                load([rootpath, '/MVPA_', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_radial.mat'])
    
                %% get weights 
                if strcmp(epoch, 'STIM')
                    % for angular modulator
                    [weights_STIM] = get_weights(subj_angular,'STIM','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_weights_STIM_angular.mat'], 'weights_STIM', '-v7.3')
                    % for radial modulator
                    [weights_STIM] = get_weights(subj_radial,'STIM','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_weights_STIM_radial.mat'], 'weights_STIM', '-v7.3')
                elseif strcmp(epoch, 'DELAY')
                    [weights_DELAY] = get_weights(subj_angular,'DELAY','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_weights_DELAY_angular.mat'], 'weights_DELAY', '-v7.3')
                    % for radial modulator
                    [weights_DELAY] = get_weights(subj_radial,'DELAY','oriRef','runs',class_args);
                    save([outpath, hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_weights_DELAY_radial.mat'], 'weights_DELAY', '-v7.3')
                end
            end
        end % of cycling over ROIs
    end % of cyclying over hemis
end % of cyclying subjects



                
% Plot reconstruct spatial map for the model prediction
% Author: Zoe Duan
% Date: 08/14/23

clear; clc;

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

% define parameters
subjIDs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
n = length(subjIDs);

% define hemisphere you want to analysis
HEMIs = {'bilat'};     % HEMIs = ['bilat', 'lh', 'rh'];
% define ROIs you want to draw
ROIs = {'V1'};

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ploting spatial reconstruction maps for modulated gratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        load([datapath, 'recon_angular.mat']);
        load([datapath, 'recon_radial.mat']);

        if sii == 1
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
    save([outpath, 'recon_angular.mat'], 'recon_angular_allSubjs', '-v7.3');
    save([outpath, 'recon_radial.mat'], 'recon_radial_allSubjs', '-v7.3');
end


%% plot
% define modulators
moduls = {'radial', 'angular'};

% define model output source
model_source = 'maxLev';
% model_source = 'avg';

% define figure formats
% fig_form = 'separate';       % separate each orientation condition
fig_form = 'combine';       % rotate all orientation conditions to align them centered at zero

% define line-fitting form
% fit_form = '';
% fit_form = 'fitline_all';   % fit lines based on all pixels on the image
fit_form = 'fitline_stim';  % fit lines only based on pixels around the stimulus

% define paths of model output
datapath = [mainpath, exptSystem];
if strcmp(model_source, 'maxLev')
    newpath = [datapath, '/results_spaceRecon_model_maxLev/'];
elseif strcmp(model_source, 'avg')
    newpath = [datapath, '/results_spaceRecon_model_avg/'];
end
% define paths to save figures
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

% load recon 
load([datapath, 'recon_angular.mat']);
load([datapath, 'recon_radial.mat']);


% Cycle over hemispheres
for hii = 1:length(HEMIs)
    hemi = HEMIs{hii};

    % Cycle over ROIs
    for rii = 1:length(ROIs)
        roi_name = ROIs{rii};

        % create recon maps for current ROI
        recon_curr = cell(1, length(moduls), n);

        % plot spatial map for current ROI
        if n == 1
            % delete empty [] in recon, because lack of ROI
            emptyCells = all(cellfun(@isempty, recon_angular), 1);
            recon_angular = recon_angular(:, ~emptyCells);
            emptyCells = all(cellfun(@isempty, recon_radial), 1);
            recon_radial = recon_radial(:, ~emptyCells);

            % create recon maps
            recon_curr(1,1,:) = recon_radial(hii,rii,:);
            recon_curr(1,2,:) = recon_angular(hii,rii,:);
            
            % save subject-wised plots
            final_fig = [outpath_fig, 'subj', num2str(subjIDs(1)), '/'];
            if ~exist(final_fig, 'dir')
                mkdir(final_fig)
            end
            savePath = [final_fig, hemi, '_', roi_name, '_spatialMaps_model_subj', num2str(subjIDs(1))];
        else         
            recon_curr(1,1,:) = recon_radial_allSubjs(hii,rii,:);
            recon_curr(1,2,:) = recon_angular_allSubjs(hii,rii,:);

            % save allSubjs plots
            savePath = [outpath_fig, hemi, '_', roi_name, '_spatialMaps_model'];
        end

        if strcmp(fig_form, 'separate')
            plot_recon_sep(recon_curr, para, [], rii, roi_name, savePath, fit_form);
        elseif strcmp(fig_form, 'combine')
            plot_recon_com(recon_curr, para, [], rii, roi_name, savePath, fit_form);
        end

    end % of cyclying ROIs
end % of cyclying hemis





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot spatial reconstruction maps for lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine all subjects data together if n>1
if n > 1
    % initializa variables
    recon_lines_allSubjs = [];

    % Cycle over subj
    for sii = 1:n
        subjID = subjIDs(sii);

        % define path for loading data
        datapath = [newpath, 'subj', num2str(subjID), '/'];

        % load recon for weights from classifier
        load([datapath, 'recon_lines.mat']);

        if sii == 1
            recon_lines_allSubjs = recon_lines;
        else
            recon_lines_allSubjs = cat(3, recon_lines_allSubjs, recon_lines);
        end
    end % of cyclying over subj

    % save combined results for all subjects
    outpath = [newpath, 'allSubjs/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % save recon results into disk
    save([outpath, 'recon_lines.mat'], 'recon_lines_allSubjs', '-v7.3');
end

%% plot
% define model output source
model_source = 'maxLev';
% model_source = 'avg';

% define figure formats
% fig_form = 'separate';       % separate each orientation condition
fig_form = 'combine';       % rotate all orientation conditions to align them centered at zero

% define line-fitting form
% fit_form = '';
% fit_form = 'fitline_all';   % fit lines based on all pixels on the image
fit_form = 'fitline_stim';  % fit lines only based on pixels around the stimulus

% define paths of model output
datapath = [mainpath, exptSystem];
if strcmp(model_source, 'maxLev')
    newpath = [datapath, '/results_spaceRecon_model_maxLev/'];
elseif strcmp(model_source, 'avg')
    newpath = [datapath, '/results_spaceRecon_model_avg/'];
end

% define paths to save figures
if strcmp(fig_form, 'separate')
    para.mapSize = [0 0 1200 600];
    outpath_fig = [newpath, 'Figures_separate/'];
elseif strcmp(fig_form, 'combine')
    para.mapSize = [0 0 600 600];
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

% load recon 
load([datapath, 'recon_lines.mat']);

% Cycle over hemispheres
for hii = 1:length(HEMIs)
    hemi = HEMIs{hii};

    % Cycle over ROIs
    for rii = 1:length(ROIs)
        roi_name = ROIs{rii};

        % create recon maps for current ROI
        recon_curr = cell(1, 1, n);

        % plot spatial map for current ROI
        if n == 1
            % delete empty [] in recon, because lack of ROI
            emptyCells = all(cellfun(@isempty, recon_lines), 1);
            recon_lines = recon_lines(:, ~emptyCells);

            % create recon maps
            recon_curr(1,1,:) = recon_lines(hii,rii,:);

            % save subject-wised plots
            final_fig = [outpath_fig, 'subj', num2str(subjIDs(1)), '/'];
            if ~exist(final_fig, 'dir')
                mkdir(final_fig)
            end
            savePath = [final_fig, hemi, '_', roi_name, '_spatialMaps_lines_subj', num2str(subjIDs(1))];
        else
            recon_curr(1,1,:) = recon_lines_allSubjs(hii,rii,:);

            % save allSubjs plots
            savePath = [outpath_fig, hemi, '_', roi_name, '_spatialMaps_lines'];
        end

        if strcmp(fig_form, 'separate')
            plot_recon_sep(recon_curr, para, [], rii, roi_name, savePath, fit_form);
        elseif strcmp(fig_form, 'combine')
            plot_recon_com(recon_curr, para, [], rii, roi_name, savePath, fit_form);
        end

    end % of cyclying ROIs
end % of cyclying hemis






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot spatial reconstruction maps for one dot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine all subjects data together if n>1
if n > 1
    % initializa variables
    recon_dots_allSubjs = [];

    % Cycle over subj
    for sii = 1:n
        subjID = subjIDs(sii);

        % define path for loading data
        datapath = [newpath, 'subj', num2str(subjID), '/'];

        % load recon for weights from classifier
        load([datapath, 'recon_1dot.mat']);

        if sii == 1
            recon_dots_allSubjs = recon_dots;
        else
            recon_dots_allSubjs = cat(3, recon_dots_allSubjs, recon_dots);
        end
    end % of cyclying over subj

    % save combined results for all subjects
    outpath = [newpath, 'allSubjs/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    % save recon results into disk
    save([outpath, 'recon_1dot.mat'], 'recon_dots_allSubjs', '-v7.3');
end


%% plot
% define model output source
model_source = 'maxLev';
% model_source = 'avg';

% define figure formats
% fig_form = 'separate';       % separate each orientation condition
fig_form = 'combine';       % rotate all orientation conditions to align them centered at zero

% define dot-fitting form
% fit_form = '';
% fit_form = 'fitline_all';   % fit lines based on all pixels on the image
fit_form = 'fitline_stim';  % fit lines only based on pixels around the stimulus

% define paths of model output
datapath = [mainpath, exptSystem];
if strcmp(model_source, 'maxLev')
    newpath = [datapath, '/results_spaceRecon_model_maxLev/'];
elseif strcmp(model_source, 'avg')
    newpath = [datapath, '/results_spaceRecon_model_avg/'];
end

% define paths to save figures
if strcmp(fig_form, 'separate')
    para.mapSize = [0 0 1200 600];
    outpath_fig = [newpath, 'Figures_separate/'];
elseif strcmp(fig_form, 'combine')
    para.mapSize = [0 0 600 600];
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

% load recon 
load([datapath, 'recon_1dot.mat']);

% Cycle over hemispheres
for hii = 1:length(HEMIs)
    hemi = HEMIs{hii};

    % Cycle over ROIs
    for rii = 1:length(ROIs)
        roi_name = ROIs{rii};

        % create recon maps for current ROI
        recon_curr = cell(1, 1, n);

        % plot spatial map for current ROI
        if n == 1
            % delete empty [] in recon, because lack of ROI
            emptyCells = all(cellfun(@isempty, recon_dots), 1);
            recon_dots = recon_dots(:, ~emptyCells);

            % create recon maps
            recon_curr(1,1,:) = recon_dots(hii,rii,:);

            % save subject-wised plots
            final_fig = [outpath_fig, 'subj', num2str(subjIDs(1)), '/'];
            if ~exist(final_fig, 'dir')
                mkdir(final_fig)
            end
            savePath = [final_fig, hemi, '_', roi_name, '_spatialMaps_1dot_subj', num2str(subjIDs(1))];
        else
            recon_curr(1,1,:) = recon_dots_allSubjs(hii,rii,:);

            % save allSubjs plots
            savePath = [outpath_fig, hemi, '_', roi_name, '_spatialMaps_1dot'];
        end

        if strcmp(fig_form, 'separate')
            plot_recon_sep(recon_curr, para, [], rii, roi_name, savePath, fit_form);
        elseif strcmp(fig_form, 'combine')
            plot_recon_com(recon_curr, para, [], rii, roi_name, savePath, fit_form);
        end

    end % of cyclying ROIs
end % of cyclying hemis









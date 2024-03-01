% plot MVPA decoding results for WMformat project 
% for the main WM task, control task, and cross-task results
% written by Zoe Duan 
% Date: 09/08/23

clear; clc; 

%% define parameters and paths
% define expname
% expname = 'WM';
% expname = 'control';
expname = 'cross';

% define modulators
moduls = {'angular', 'radial'};
% define decoding type
decodingType = {'Within', 'Cross'};

% define subj information
subjIDs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
n = length(subjIDs);

% define paths
mainpath = input('Input the path for the parent folder', 's');
% mainpath = '/System/Volumes/Data/d/DATC/datc/formatwm_ZD/public/';

if strcmp(expname, 'WM') 
    datapath_fmri = [mainpath, 'results_decoding_WMTask/'];
    outpath_fig = [mainpath, 'results_decoding_WMTask/Figures/'];
    % define epoch
    epoch = 'DELAY';
elseif strcmp(expname, 'control')
    datapath_fmri = [mainpath, 'results_decoding_controlTask/'];
    outpath_fig = [mainpath, 'results_decoding_controlTask/Figures/'];
    % define epoch
    epoch = 'STIM';
elseif strcmp(expname, 'cross')
    datapath_fmri = [mainpath, 'results_decoding_crossTask/'];
    outpath_fig = [mainpath, 'results_decoding_crossTask/Figures/'];
    % define epoch
    epoch = 'DELAY';
end
if ~exist(outpath_fig, 'dir')
   mkdir(outpath_fig)
end


% define hemisphere you want to analysis
HEMIs = {'bilat'};     % HEMIs = ['bilat', 'lh', 'rh'];
% define ROIs you want to draw
% ROIs = {'V1', 'V2', 'V3', 'V3AB', 'IPS0_IPS1'};
ROIs = {'V1', 'V2', 'V3', 'V3AB', 'IPS0_IPS1', 'IPS2_IPS3', 'sPCS', 'iPCS'};

% define parameters
condcolors = [190 0 110; 0 110 190]./255;
if length(ROIs) == 8
    fig_size = [0 0 1800 1000];
else
    fig_size = [0 0 1800 400];
end


%% combine all subject results together
if n > 1
    % Cycle over hemispheres
    for hii = 1:length(HEMIs)
        hemi = HEMIs{hii};

        % Cycle over ROIs
        for rii = 1:length(ROIs)
            roi_name = ROIs{rii};

            % initialize variables
            results_angular_allSubj = [];
            results_a2r_allSubj = [];
            results_radial_allSubj = [];
            results_r2a_allSubj = [];

            % Cycle over subj
            for sii = 1:n
                subjID = subjIDs(sii);

                %% for within angular modulator
                filename = [datapath_fmri, 'subj', num2str(subjID), '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_', epoch, '_angular.mat'];
                if ~exist(filename, 'file')
                    sprintf([hemi, '_', roi_name, '_', 'subj', num2str(subjID), 'results_within_angular does not exist'])
                else
                    % load data for current subj
                    load(filename);
                    if strcmp(epoch, 'STIM')
                        results_angular = results_STIM_angular;
                    elseif strcmp(epoch, 'DELAY')
                        results_angular = results_DELAY_angular;
                    end
                    % get the total performance
                    results_angular_currSubj = results_angular.total_perf;
                    % add up all subjs data
                    results_angular_allSubj = [results_angular_allSubj, results_angular_currSubj];
                end

                %% for cross a2r modulator
                filename = [datapath_fmri, 'subj', num2str(subjID), '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_', epoch, '_a2r.mat'];
                if ~exist(filename, 'file')
                    sprintf([hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_a2r.mat does not exist'])
                else
                    % load data for current subj
                    load(filename);
                    if strcmp(epoch, 'STIM')
                        results_a2r = results_STIM_a2r;
                    elseif strcmp(epoch, 'DELAY')
                        results_a2r = results_DELAY_a2r;
                    end
                    % get the total performance
                    results_a2r_currSubj = results_a2r.total_perf;
                    % add up all subjs data
                    results_a2r_allSubj = [results_a2r_allSubj, results_a2r_currSubj];
                end

                %% for within radial modulator
                filename = [datapath_fmri, 'subj', num2str(subjID), '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_', epoch, '_radial.mat'];
                if ~exist(filename, 'file')
                    sprintf([hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_within_radial.mat does not exist'])
                else
                    % load data for current subj
                    load(filename);
                    if strcmp(epoch, 'STIM')
                        results_radial = results_STIM_radial;
                    elseif strcmp(epoch, 'DELAY')
                        results_radial = results_DELAY_radial;
                    end
                    % get the total performance
                    results_radial_currSubj = results_radial.total_perf;
                    % add up all subjs data
                    results_radial_allSubj = [results_radial_allSubj, results_radial_currSubj];
                end

                %% for cross r2a modulator
                filename = [datapath_fmri, 'subj', num2str(subjID), '/', hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_', epoch, '_r2a.mat'];
                if ~exist(filename, 'file')
                    sprintf([hemi, '_', roi_name, '_', 'subj', num2str(subjID), '_results_cross_r2a.mat does not exist'])
                else
                    % load data for current subj
                    load(filename);
                    if strcmp(epoch, 'STIM')
                        results_r2a = results_STIM_r2a;
                    elseif strcmp(epoch, 'DELAY')
                        results_r2a = results_DELAY_r2a;
                    end
                    % get the total performance
                    results_r2a_currSubj = results_r2a.total_perf;
                    % add up all subjs data
                    results_r2a_allSubj = [results_r2a_allSubj, results_r2a_currSubj];
                end
            end % of cycling over subj
            
            % save results
            newpath = [datapath_fmri, 'allSubjs/'];
            if ~exist(newpath, 'dir')
                mkdir(newpath)
            end
            
            save([newpath, hemi, '_', roi_name, '_results_within_', epoch, '_angular.mat'], 'results_angular_allSubj', '-v7.3');
            save([newpath, hemi, '_', roi_name, '_results_within_', epoch, '_radial.mat'], 'results_radial_allSubj', '-v7.3');
            save([newpath, hemi, '_', roi_name, '_results_cross_', epoch, '_a2r.mat'], 'results_a2r_allSubj', '-v7.3');
            save([newpath, hemi, '_', roi_name, '_results_cross_', epoch, '_r2a.mat'], 'results_r2a_allSubj', '-v7.3');
        end % of cycling over ROIs
    end % of cycling over hemi
end


%% plot results
% Cycle over hemispheres
for hii = 1:length(HEMIs)
    hemi = HEMIs{hii};

    % Cycle over ROIs
    for rii = 1:length(ROIs)
        roi_name = ROIs{rii};

        %% load data
        % define datapath for loading data 
        if n == 1
            datapath = [datapath_fmri, 'subj', num2str(subjIDs(1)), '/'];
        else
            datapath = [datapath_fmri, 'allSubjs/'];
        end
        % load data
        load([datapath, hemi, '_', roi_name, '_results_within_', epoch, '_angular.mat']);
        load([datapath, hemi, '_', roi_name, '_results_within_', epoch, '_radial.mat']);
        load([datapath, hemi, '_', roi_name, '_results_cross_', epoch, '_a2r.mat']);
        load([datapath, hemi, '_', roi_name, '_results_cross_', epoch, '_r2a.mat']);
        
        if n == 1
            results_angular_allSubj = results_angular;
            results_radial_allSubj = results_radial;
            results_a2r_allSubj = results_a2r;
            results_r2a_allSubj = results_r2a;
        end


        %% calculate mean and SEM for decoding results 
        % calculate mean and std cross subj
        results_angular_mean = mean(results_angular_allSubj);
        results_angular_SEM = std(results_angular_allSubj) ./ sqrt(n);
        results_a2r_mean = mean(results_a2r_allSubj);
        results_a2r_SEM = std(results_a2r_allSubj) ./ sqrt(n);

        results_radial_mean = mean(results_radial_allSubj);
        results_radial_SEM = std(results_radial_allSubj) ./ sqrt(n);
        results_r2a_mean = mean(results_r2a_allSubj);
        results_r2a_SEM = std(results_r2a_allSubj) ./ sqrt(n);


        %% plot all conditions on the same picture
        % combine different condtions together
        mean_angular = [results_angular_mean; results_a2r_mean]; 
        mean_radial = [results_radial_mean; results_r2a_mean];
        mean_all = [mean_radial, mean_angular];
        SEM_angular = [results_angular_SEM; results_a2r_SEM];
        SEM_radial = [results_radial_SEM; results_r2a_SEM];
        SEM_all = [SEM_radial, SEM_angular];

        figure(1)
        if length(ROIs)== 8
            subplot(2, length(ROIs)/2, rii)
        else
            subplot(1, length(ROIs), rii)
        end
        fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 20;
        set(gcf,'Position', fig_size)
        set(gcf, 'DefaultAxesFontName', 'Arial');

        hold on

        % plot bar
        b = bar(mean_all, 'grouped', 'BarWidth', 1);
        % set bar color
        set(b(1), 'FaceColor', condcolors(1, :));
        set(b(2), 'FaceColor', condcolors(2, :));
        % set bar distance between groups
        n_groups = length(moduls);
        n_barsPerGroup = length(decodingType);
        
        % Calculate the width of each group
        
        groupWidth = min(0.8, n_groups/(n_groups+1.5));
        % Loop through each modul to add error bars
        for i = 1:n_groups
            x = (1:n_barsPerGroup) - groupWidth/2 + (2*i-1) * groupWidth / (2*n_groups);
            errorbar(x, mean_all(:,i), SEM_all(:,i), 'k', 'linestyle', 'none', 'linewidth', 2);
        end

        % plot individual datapoints
        scatter(repmat(b(1).XEndPoints(1), length(results_radial_allSubj), 1),results_radial_allSubj, ...
            60,'MarkerFaceColor','none','MarkerEdgeColor',condcolors(1, :),'MarkerEdgeAlpha',.5,'LineWidth',1,'jitter','on', 'jitterAmount',0.05);
        scatter(repmat(b(2).XEndPoints(1), length(results_angular_allSubj), 1),results_angular_allSubj, ...
            60,'MarkerFaceColor','none','MarkerEdgeColor',condcolors(2, :),'MarkerEdgeAlpha',.5,'LineWidth',1,'jitter','on', 'jitterAmount',0.05);
        scatter(repmat(b(1).XEndPoints(2), length(results_r2a_allSubj), 1),results_r2a_allSubj, ...
            60,'MarkerFaceColor','none','MarkerEdgeColor',condcolors(1, :),'MarkerEdgeAlpha',.5,'LineWidth',1,'jitter','on', 'jitterAmount',0.05);
        scatter(repmat(b(2).XEndPoints(2), length(results_a2r_allSubj), 1),results_a2r_allSubj, ...
            60,'MarkerFaceColor','none','MarkerEdgeColor',condcolors(2, :),'MarkerEdgeAlpha',.5,'LineWidth',1,'jitter','on', 'jitterAmount',0.05);

        % draw chance level line
        plot(xlim, [1/3 1/3],'k--','LineWidth',1.25)
        
        % define xticks
        xticks(1:1:n_barsPerGroup)
        xticklabels(decodingType);

        % set others
        % fix the issue of presenting IPS0_IPS1 and IPS2_IPS3
        if strcmp(roi_name, 'IPS0_IPS1')
            roi_name = 'IPS0/1';
        elseif strcmp(roi_name, 'IPS2_IPS3')
            roi_name = 'IPS2/3';
        end
%         title(roi_name)
        if rii == 1
            ylabel('Decoding accuracy')
        end
        ylim([0.0, 0.8])
        yticks(0.1*(0:2:9))

        legend('off');    % show legend for the last plot only

        if rii == length(ROIs)
            % show legend
%             legend('angular modulator', 'radial modulator');
%             legend('boxoff')
%             legend('Location','best')
            % last ROI has been plotted
            if n == 1
                % save subject-wise results
                saveas(fig, [outpath_fig, 'decoding_results_', epoch, '_', expname, 'Task_subj', num2str(subjIDs(1)), '.png']);
                set(fig, 'Renderer', 'painters'); % make text ediable
                saveas(fig, [outpath_fig, 'decoding_results_', epoch, '_', expname, 'Task_subj', num2str(subjIDs(1)), '.eps'], 'epsc');
            else
                saveas(fig, [outpath_fig, 'decoding_results_', epoch, '_', expname, 'Task_', num2str(n), 'subjs_', num2str(length(ROIs)),'ROIs.png']);
                set(fig, 'Renderer', 'painters'); % make text ediable
                saveas(fig, [outpath_fig, 'decoding_results_', epoch, '_', expname, 'Task_', num2str(n), 'subjs_', num2str(length(ROIs)),'ROIs.eps'], 'epsc');
            end
        end
        pause(1)

    end % of cycling over ROIs
end % of cycling over hemi



% Plot filtered responses and reconstruction fidelity for both tasks
% Author: Zoe Duan 
% Date: 08/30/23

clear; clc;

% add paths
addpath(genpath('/System/Volumes/Data/d/DATA/home/zoe/vistasoft_ts/'));
addpath(genpath('/System/Volumes/Data/d/DATA/home/zoe/preproc_mFiles/'));
addpath('/Users/zoe/Github/princeton-mvpa-toolbox'); 
mvpa_add_paths;

% define expname
expname = 'WM';
% expname = 'control';

% define line filter form
% filter_form = 'full';     % use all pixels of the recon image
filter_form = 'stim';     % use pixels within the stimulu

% define subjs info
subjIDs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
n = length(subjIDs);

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

% define path to save figures
outpath_fig = [datapath_main, 'Figures_fidelity_filterForm_', filter_form, '/'];
if ~exist(outpath_fig, 'dir')
   mkdir(outpath_fig)
end

%% combine all subjects data together if n>1
if n > 1
    % initializa variables
    filterResp_angular_allSubjs = [];
    filterResp_radial_allSubjs = [];
    fidelity_angular_allSubjs = [];
    fidelity_radial_allSubjs = [];

    % Cycle over subj
    for sii = 1:n
        subjID = subjIDs(sii);

        % define path for loading data
        datapath = [datapath_main, 'subj', num2str(subjID), '/'];
        % load data
        load([datapath, 'filterResp_angular_', epoch, '_filterForm_', filter_form, '.mat']);
        load([datapath, 'filterResp_radial_', epoch, '_filterForm_', filter_form, '.mat']);

        load([datapath, 'fidelity_angular_', epoch, '_filterForm_', filter_form, '.mat']);
        load([datapath, 'fidelity_radial_', epoch, '_filterForm_', filter_form, '.mat']);

        % combine all subj data together
        if isempty(filterResp_angular_allSubjs)
            filterResp_angular_allSubjs = filterResp_angular;
            filterResp_radial_allSubjs = filterResp_radial;

            fidelity_angular_allSubjs = fidelity_angular;
            fidelity_radial_allSubjs = fidelity_radial;
        else
            filterResp_angular_allSubjs = cat(3, filterResp_angular_allSubjs, filterResp_angular);
            filterResp_radial_allSubjs = cat(3, filterResp_radial_allSubjs, filterResp_radial);

            fidelity_angular_allSubjs = cat(3, fidelity_angular_allSubjs, fidelity_angular);
            fidelity_radial_allSubjs = cat(3, fidelity_radial_allSubjs, fidelity_radial);
        end
    end % of cyclying over subj

    % save combined results for all subjects
    outpath = [datapath_main, 'allSubjs/'];
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end

    save([outpath, 'filterResp_angular_', epoch, '_filterForm_', filter_form, '.mat'], 'filterResp_angular_allSubjs', '-v7.3');
    save([outpath, 'filterResp_radial_', epoch, '_filterForm_', filter_form, '.mat'], 'filterResp_radial_allSubjs', '-v7.3');

    save([outpath, 'fidelity_angular_', epoch, '_filterForm_', filter_form, '.mat'], 'fidelity_angular_allSubjs', '-v7.3');
    save([outpath, 'fidelity_radial_', epoch, '_filterForm_', filter_form, '.mat'], 'fidelity_radial_allSubjs', '-v7.3');

end



%% plot filtered responses
% define hemisphere you want to analysis
HEMIs = {'bilat'};     % HEMIs = ['bilat', 'lh', 'rh'];
% define ROIs you want to draw
% ROIs = {'V1', 'V2', 'V3', 'V3AB', 'IPS0_IPS1'};
ROIs = {'V1', 'V2', 'V3', 'V3AB', 'IPS0_IPS1', 'IPS2_IPS3', 'sPCS', 'iPCS'};      % put these ROIs in supplementary

% define parameters for plots
condcolors = [190 0 110; 0 110 190]./255;
if length(ROIs)==8
    fig_size_filteredResp = [0 0 2500 500];
    fig_size_fidelity = [0 0 1500 800];
else
    fig_size_filteredResp = [0 0 2200 500];
    fig_size_fidelity = [0 0 1000 500];
end

y_lim = [-2, 2];
n_oris = 180;


% define datapath for loading data
if n == 1
    datapath = [datapath_main, 'subj', num2str(subjIDs(1)), '/'];
else
    datapath = [datapath_main, 'allSubjs/'];
end
% load recon for real data
load([datapath, 'filterResp_angular_', epoch, '_filterForm_', filter_form, '.mat']);
load([datapath, 'filterResp_radial_', epoch, '_filterForm_', filter_form, '.mat']);

% Cycle over hemispheres
for hii = 1:length(HEMIs)
    hemi = HEMIs{hii};

    % Cycle over ROIs
    for rii = 1:length(ROIs)
        roi_name = ROIs{rii};
        % fix the issue of presenting IPS0_IPS1 and IPS2_IPS3
        if strcmp(roi_name, 'IPS0_IPS1')
            roi_name = 'IPS0/1';
        elseif strcmp(roi_name, 'IPS2_IPS3')
            roi_name = 'IPS2/3';
        end

        if n == 1
            filterResp_angular_mean = mean(filterResp_angular{hii,rii}, 2)';
            filterResp_angular_SEM = (std(filterResp_angular{hii,rii}, [], 2) ./ sqrt(n))';

            filterResp_radial_mean = mean(filterResp_radial{hii,rii}, 2)';
            filterResp_radial_SEM = (std(filterResp_radial{hii,rii}, [], 2) ./ sqrt(n))';
        else
            % calculate mean and std across subj
            r_angular = [cellfun(@(x) {x(:,:)}, filterResp_angular_allSubjs(hii,rii,:))];
            rcat_angular = squeeze(cat(3,r_angular{:}))';
            filterResp_angular_mean = mean(rcat_angular, 1);
            filterResp_angular_SEM = std(rcat_angular, [], 1) ./ sqrt(n);

            r_radial = [cellfun(@(x) {x(:,:)}, filterResp_radial_allSubjs(hii,rii,:))];
            rcat_radial = squeeze(cat(3,r_radial{:}))';
            filterResp_radial_mean = mean(rcat_radial, 1);
            filterResp_radial_SEM = std(rcat_radial, [], 1) ./ sqrt(n);
        end

        % plot both modulators on the same picture
        figure(1)
        f=subplot(1, length(ROIs), rii);

        hold on
        % plot mean
        plot(1:n_oris,filterResp_radial_mean(1:n_oris),'Color',condcolors(1,:),'LineWidth',2, 'HandleVisibility','off')
        plot(1:n_oris,filterResp_angular_mean(1:n_oris),'Color',condcolors(2,:),'LineWidth',2, 'HandleVisibility','off')
        % plot errorbar
        btwn_fill = [filterResp_radial_mean(1:n_oris) + filterResp_radial_SEM(1:n_oris), fliplr(filterResp_radial_mean(1:n_oris) - filterResp_radial_SEM(1:n_oris))];
        fill_xs = [1:n_oris, fliplr(1:n_oris)];
        fill(fill_xs,btwn_fill,condcolors(1,:),'linestyle','none','facealpha',0.2, 'HandleVisibility','off');
        btwn_fill = [filterResp_angular_mean(1:n_oris) + filterResp_angular_SEM(1:n_oris), fliplr(filterResp_angular_mean(1:n_oris) - filterResp_angular_SEM(1:n_oris))];
        fill_xs = [1:n_oris, fliplr(1:n_oris)];
        fill(fill_xs,btwn_fill,condcolors(2,:),'linestyle','none','facealpha',0.2, 'HandleVisibility','off');

        % define xticks
        xticks(0:n_oris/2:n_oris)
        xticklabels(-90:90:90)
        % set ylim
        ylim(y_lim);

        % draw center line
        plot([n_oris/2 n_oris/2],ylim,'k','LineWidth',1.25,'HandleVisibility','off')

        % set others
        fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 20;
        set(gcf,'Position',fig_size_filteredResp)
        %             f.Title.String = roi_name;
        %             f.Title.FontSize = 14;
        if rii == 1
            ylabel('Filtered responses');
        end

        xlabel('Orientation filters')
        legend('off')
        legend('boxoff')
        %             legend('angular modulator', 'radial modulator');
        %             legend('Location','best')

        if rii == length(ROIs)
            % add subtitle for the figure
            %                 suptitle('Filtered responses for DELAY period', 'FontSize', 20, 'FontWeight', 'bold');
            % save figure when last ROI has been plotted
            if n == 1
                % save subject-wised plots
                final_fig = [outpath_fig, 'subj', num2str(subjIDs(1)), '/'];
                if ~exist(final_fig, 'dir')
                    mkdir(final_fig)
                end
                savePath = [final_fig, 'filterResp', epoch, '_real_subj', num2str(subjIDs(1)), '_', num2str(length(ROIs)), 'ROIs'];
            else
                % save allSubjs plots
                savePath = [outpath_fig, 'filterResp', epoch, '_real_', num2str(length(ROIs)), 'ROIs'];
            end
            saveas(fig, [savePath, '.png']);
            set(fig, 'Renderer', 'painters'); % make text ediable
            saveas(fig, [savePath, '.eps'], 'epsc');
        end
        pause(1)

    end % of cycling over ROIs
end % of cycling over hemi




%% plot reconstruction fidelity
% define datapath for loading data
if n == 1
    datapath = [datapath_main, 'subj', num2str(subjIDs(1)), '/'];
else
    datapath = [datapath_main, 'allSubjs/'];
end
% load fidelity for all subjects
load([datapath, 'fidelity_angular_', epoch, '_filterForm_', filter_form, '.mat']);
load([datapath, 'fidelity_radial_', epoch, '_filterForm_', filter_form, '.mat']);

% initialize variables
mean_angular = [];
SEM_angular = [];
mean_radial = [];
SEM_radial = [];

% Cycle over hemispheres
for hii = 1:length(HEMIs)
    hemi = HEMIs{hii};

    % Cycle over ROIs
    for rii = 1:length(ROIs)
        roi_name = ROIs{rii};
        % fix the issue of presenting IPS0_IPS1 and IPS2_IPS3
        if strcmp(roi_name, 'IPS0_IPS1')
            roi_name = 'IPS0/1';
        elseif strcmp(roi_name, 'IPS2_IPS3')
            roi_name = 'IPS2/3';
        end
        ROIs{rii} = roi_name;

        if n == 1
            fidelity_angular_mean = mean(fidelity_angular{hii,rii});
            fidelity_angular_SEM = std(fidelity_angular{hii,rii}) ./ sqrt(n);

            fidelity_radial_mean = mean(fidelity_radial{hii,rii});
            fidelity_radial_SEM = std(fidelity_radial{hii,rii}) ./ sqrt(n);
        else
            % calculate mean and std across subj
            rcat_angular = squeeze(cat(3,fidelity_angular_allSubjs{hii,rii,:}));
            fidelity_angular_mean = mean(rcat_angular);
            fidelity_angular_SEM = std(rcat_angular) ./ sqrt(n);

            rcat_radial = squeeze(cat(3,fidelity_radial_allSubjs{hii,rii,:}));
            fidelity_radial_mean = mean(rcat_radial);
            fidelity_radial_SEM = std(rcat_radial) ./ sqrt(n);
        end

        % combine results across ROIs together
        if isempty(mean_angular)
            mean_angular = fidelity_angular_mean;
            SEM_angular = fidelity_angular_SEM;
            mean_radial = fidelity_radial_mean;
            SEM_radial = fidelity_radial_SEM;
        else
            mean_angular = [mean_angular; fidelity_angular_mean];
            SEM_angular = [SEM_angular; fidelity_angular_SEM];
            mean_radial = [mean_radial; fidelity_radial_mean];
            SEM_radial = [SEM_radial; fidelity_radial_SEM];
        end
        % group angular and radial results together
        mean_final = [mean_radial, mean_angular];
        SEM_final = [SEM_radial, SEM_angular];

    end % of cycling over ROIs
end % of cycling over hemi

% plot both modulators on the same picture
figure;
hold on

% plot bar
b = bar(mean_final, 'grouped', 'BarWidth', 1);
set(b(1), 'FaceColor', condcolors(1, :));
set(b(2), 'FaceColor', condcolors(2, :));

% plot error bar
% Calculate the number of ROIs and groups
n_ROIs = size(mean_final, 1);
n_moduls = size(mean_final, 2);

% Calculate the width of each group
groupWidth = min(0.8, n_moduls/(n_moduls+1.5));

% Loop through each modul to add error bars
for i = 1:n_moduls
    x = (1:n_ROIs) - groupWidth/2 + (2*i-1) * groupWidth / (2*n_moduls);
    errorbar(x, mean_final(:,i), SEM_final(:,i), 'k', 'linestyle', 'none');
end

% define xticks
xticks(1:n_ROIs)
xticklabels(ROIs)

% set others
fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 20;
set(gcf,'Position',fig_size_fidelity)
ylabel('Fidelity (a.u.)');


%     title('Reconstruction fidelity for DELAY period', 'FontSize', 20, 'FontWeight', 'bold');
legend('off')
%     legend('angular modulator', 'radial modulator');
%     legend('boxoff')
%     legend('Location','best')

% save figure
if n == 1
    % save subject-wised plots
    final_fig = [outpath_fig, 'subj', num2str(subjIDs(1)), '/'];
    if ~exist(final_fig, 'dir')
        mkdir(final_fig)
    end
    savePath = [final_fig, 'fidelity_', epoch, '_real_subj', num2str(subjIDs(1)), '_', num2str(length(ROIs)), 'ROIs'];
else
    % save allSubjs plots
    savePath = [outpath_fig, 'fidelity_', epoch, '_real_', num2str(length(ROIs)), 'ROIs'];
end
saveas(fig, [savePath, '.png']);
set(fig, 'Renderer', 'painters'); % make text ediable
saveas(fig, [savePath, '.eps'], 'epsc');






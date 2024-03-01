% generate model output for angular and radial stimuli
% You need to download the following matlab tools and addpaths
% mrTools: https://github.com/justingardner/mrTools
% mgl: https://github.com/justingardner/mgl
% matlabPyrTools: https://github.com/LabForComputationalVision/matlabPyrTools
% stimulusVignetting: https://github.com/elimerriam/stimulusVignetting 
% Author: Zoe Duan
% Date: 05/03/2023

% add paths 
addpath(genpath('/Users/zoe/Github/mrTools/'));
addpath(genpath('/Users/zoe/Github/mgl/'));
addpath(genpath('/Users/zoe/Github/matlabPyrTools/'));
addpath(genpath('/Users/zoe/Github/stimulusVignetting/'));


clear; clc

% define the system
exptSystem = 'scanRoom';

% define paths
mainpath = input('Input the path for the parent folder', 's');
% mainpath = '/System/Volumes/Data/d/DATC/datc/formatwm_ZD/public/';
datapath = [mainpath, 'models/'];
outpath = [mainpath, 'modelResults/'];
if ~exist(outpath, 'dir')
   mkdir(outpath)
end

% load stimuli saved by create_stimScreen.m
% whose shape is (x, y, nTrig, nPhase, nOri)
load([datapath, 'final_stimScreen_angular_', exptSystem, '.mat']);
load([datapath, 'final_stimScreen_radial_', exptSystem, '.mat']);

% define oris you want to simulate from the model
% oris = [0, 90, 15, 75, 135];

% define condtions
stimP.nSquare = 1;                    % # of squares; 1:square, 2:nonSquare
stimP.nTrig = 2;                      % # of trig; 1:cos, 2:sin
stimP.nPhase = 2;                     % # of phases within each trig
stimP.phase = linspace(0,2*pi,stimP.nPhase+1);     % equi-distant phases
stimP.phase(end) = [];
stimP.sfStatic = 1;                   % spatial frequency for carrier (cycles per degree)
stimP.angFreq = 5;                    % spatial frequency for modulator

% define the model parameters
% define the # of orientation bands
numOrientations = 6; 
% compute the orientation bandwidth
% oriBandwidth = 360/numOrientations;   
% define the spatial frequency bandwidth
bandwidth = 0.5;        


%% build pryamid for angular modulator
% find the dimension of image
dims=size(final_stimScreen_angular, [1,2]);
% compute the # of levels/scales based on the bandwidth and image size
numLevels = maxLevel(dims,bandwidth);
% construct quad frequency filters
[freqRespsImag,freqRespsReal,modulDiff]= makeQuadFRs(dims,numLevels,numOrientations,bandwidth);

% initialize variables
sumBands_angular = zeros([dims,stimP.nTrig,stimP.nPhase,numel(oris)]);

% Cycle over oris
for i = 1:length(oris)
    % Cycle over trigs
    for trig = 1:stimP.nTrig
        % Cycle over phases
        for phase = 1:stimP.nPhase
            % find the image
            im = final_stimScreen_angular(:,:,trig,phase,i);
            [pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
            % Cycle over spatial frequency bands
            for lev = 1:numLevels
                sumBandsLev_angular{lev}(:,:,trig,phase,i) = zeros(dims);
                % Cycle over orientation bands
                for orientation = 1:numOrientations  
                    % extract frequency response
                    thisBand = accessSteerBand(pyr,pind,numOrientations,lev,orientation);
                    bandIm = abs(thisBand).^2;
                    % sum across orientation channels
                    sumBands_angular(:,:,trig,phase,i) = squeeze(sumBands_angular(:,:,trig,phase,i))+bandIm;
                    sumBandsLev_angular{lev}(:,:,trig,phase,i) = squeeze(sumBandsLev_angular{lev}(:,:,trig,phase,i))+bandIm;
                end % of over orientation bands
            end % of cycling over levels
        end % of cycling over phases
    end % of cycling over trigs
    
    % average the response across trigs (sin/cos)
    for phase = 1:stimP.nPhase
        trigAvgBands_angular(:,:,phase,i) = squeeze(mean(sumBands_angular(:,:,:,phase,i), 3));
        for lev = 1:numLevels
            trigAvgBandsLev_angular{lev}(:,:,phase,i) = squeeze(mean(sumBandsLev_angular{lev}(:,:,:,phase,i), 3));
        end
    end
    
    % average across phases
    for trig = 1:stimP.nTrig
        phaseAvgBands_angular(:,:,trig,i) = squeeze(mean(sumBands_angular(:,:,trig,:,i), 4));
        for lev=1:numLevels
            phaseAvgBandsLev_angular{lev}(:,:,trig,i) = squeeze(mean(sumBandsLev_angular{lev}(:,:,trig,:,i), 4));
        end
    end

    % average across both trigs and phases
    avgBands_angular(:,:,i) = squeeze(mean(trigAvgBands_angular(:,:,:,i), 3));
    for lev=1:numLevels
        avgBandsLev_angular{lev}(:,:,i) = squeeze(mean(trigAvgBandsLev_angular{lev}(:,:,:,i), 3));
    end
end  % of cycling over oris

% save results
save([outpath, '/modelOutput_angular.mat'],'sumBands_angular','sumBandsLev_angular','numOrientations','bandwidth','dims','numLevels','oris','-v7.3');
save([outpath, '/avgModelOutput_angular.mat'], 'avgBands_angular', 'avgBandsLev_angular', 'numOrientations','bandwidth','dims','numLevels','oris','-v7.3');
save([outpath, '/avgModelOutput_others_angular.mat'], 'trigAvgBands_angular', 'trigAvgBandsLev_angular', 'phaseAvgBands_angular','phaseAvgBandsLev_angular', '-v7.3');


%% build pryamid for radial modulator
% find the dimension of image
dims=size(final_stimScreen_radial, [1,2]);
% compute the # of spatial frequency bands based on the bandwidth and image size
numLevels = maxLevel(dims,bandwidth);
% construct quad frequency filters
[freqRespsImag,freqRespsReal,modulDiff]= makeQuadFRs(dims,numLevels,numOrientations,bandwidth);

% initialize variables
sumBands_radial = zeros([dims,stimP.nTrig,stimP.nPhase,numel(oris)]);

% Cycle over oris
for i = 1:length(oris)
    % Cycle over trigs
    for trig = 1:stimP.nTrig
        % Cycle over phases
        for phase = 1:stimP.nPhase
            % find the image
            im = final_stimScreen_radial(:,:,trig,phase,i);
            [pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
            % Cycle over spatial frequency bands
            for lev = 1:numLevels
                sumBandsLev_radial{lev}(:,:,trig,phase,i) = zeros(dims);
                % Cycle over orientation bands
                for orientation = 1:numOrientations  
                    % extract frequency response
                    thisBand = accessSteerBand(pyr,pind,numOrientations,lev,orientation);
                    bandIm = abs(thisBand).^2;
                    % sum across orientation channels
                    sumBands_radial(:,:,trig,phase,i) = squeeze(sumBands_radial(:,:,trig,phase,i))+bandIm;
                    sumBandsLev_radial{lev}(:,:,trig,phase,i) = squeeze(sumBandsLev_radial{lev}(:,:,trig,phase,i))+bandIm;
                end % of over orientation bands
            end % of cycling over levels
        end % of cycling over phases
    end % of cycling over trigs
    
    % average the response across trigs (sin/cos)
    for phase = 1:stimP.nPhase
        trigAvgBands_radial(:,:,phase,i) = squeeze(mean(sumBands_radial(:,:,:,phase,i), 3));
        for lev = 1:numLevels
            trigAvgBandsLev_radial{lev}(:,:,phase,i) = squeeze(mean(sumBandsLev_radial{lev}(:,:,:,phase,i), 3));
        end
    end
    
    % average across phases
    for trig = 1:stimP.nTrig
        phaseAvgBands_radial(:,:,trig,i) = squeeze(mean(sumBands_radial(:,:,trig,:,i), 4));
        for lev=1:numLevels
            phaseAvgBandsLev_radial{lev}(:,:,trig,i) = squeeze(mean(sumBandsLev_radial{lev}(:,:,trig,:,i), 4));
        end
    end

    % average across both trigs and phases
    avgBands_radial(:,:,i) = squeeze(mean(trigAvgBands_radial(:,:,:,i), 3));
    for lev=1:numLevels
        avgBandsLev_radial{lev}(:,:,i) = squeeze(mean(trigAvgBandsLev_radial{lev}(:,:,:,i), 3));
    end
end  % of cycling over oris

% save results
save([outpath, '/modelOutput_radial.mat'],'sumBands_radial','sumBandsLev_radial','numOrientations','bandwidth','dims','numLevels','oris','-v7.3');
save([outpath, '/avgModelOutput_radial.mat'], 'avgBands_radial', 'avgBandsLev_radial', 'numOrientations','bandwidth','dims','numLevels','oris','-v7.3');
save([outpath, '/avgModelOutput_others_radial.mat'], 'trigAvgBands_radial', 'trigAvgBandsLev_radial', 'phaseAvgBands_radial','phaseAvgBandsLev_radial', '-v7.3');


%% Find level with maximal difference between modulators
% for vertical vs. horizontal orientations
for lev=1:numLevels
    % for angular modulator
    ver_angular = squeeze(avgBandsLev_angular{lev}(:,:,1));
    hor_angular = squeeze(avgBandsLev_angular{lev}(:,:,2));
    oriDiff_angular{lev} = ver_angular-hor_angular;

    % for radial modulator
    ver_radial = squeeze(avgBandsLev_radial{lev}(:,:,1));
    hor_radial = squeeze(avgBandsLev_radial{lev}(:,:,2));
    oriDiff_radial{lev} = ver_radial-hor_radial;

    modulDiff{lev} = oriDiff_angular{lev}-oriDiff_radial{lev};
    maxModulDiff(lev) = max(modulDiff{lev}(:));
end
[tmp, maxLev] = max(maxModulDiff);

% save results
save([outpath, '/modulDiff_maxLev.mat'], 'maxLev', 'oriDiff_angular', 'oriDiff_radial', 'modulDiff', 'maxModulDiff');



%% plot

%% plot model output results for angular modulator
% Cycle over oris
for i = 1:length(oris)
    % define the output path
    output_dir = [outpath, '/modelOutput_plots_angular_', num2str(oris(i)), 'deg/'];
    if ~exist(output_dir, 'dir')
       mkdir(output_dir)
    end

    % plot averaged results
    figure('Position', [0, 0, 960, 540]);
    fig = imagesc(avgBands_angular(:,:,i));
    colormap gray
    axis off
    saveas(fig, [output_dir, 'allLev_angular.png'])
    
    % plot results for each level
    for lev=1:numLevels
        figure('Position', [0, 0, 960, 540]);
        fig = imagesc(avgBandsLev_angular{lev}(:,:,i));
        colormap gray
        axis off
        saveas(fig, [output_dir, 'lev', num2str(lev),'_angular.png'])
    end
end % of cycling over 
close all


%% plot model output results for radial modulator
for i = 1:length(oris)
    % define the output path
    output_dir = [outpath, '/modelOutput_plots_radial_', num2str(oris(i)), 'deg/'];
    if ~exist(output_dir, 'dir')
       mkdir(output_dir)
    end

    % plot averaged results
    figure('Position', [0, 0, 960, 540]);
    fig = imagesc(avgBands_radial(:,:,i));
    colormap gray
    axis off
    saveas(fig, [output_dir, 'allLev_radial.png'])
    
    % plot results for each level
    for lev=1:numLevels
        figure('Position', [0, 0, 960, 540]);
        fig = imagesc(avgBandsLev_radial{lev}(:,:,i));
        colormap gray
        axis off
        saveas(fig, [output_dir, 'lev', num2str(lev),'_radial.png'])
    end
end % of cycling over 
close all



%% plot the response to vertical and horizontal stimuli for each trig at the maximum level
% for angular modulator
for trig = 1:1:stimP.nTrig
    % for vertical
    figure('Position', [0, 0, 960, 540]);
    fig = imagesc(phaseAvgBandsLev_angular{maxLev}(:,:,trig,1));
    colormap gray
    axis off
    saveas(fig, [outpath, ['/Vertical_trig', num2str(trig), '_maxLev_angular.png']])

    figure('Position', [0, 0, 960, 540]);
    % for horizontal
    fig = imagesc(phaseAvgBandsLev_angular{maxLev}(:,:,trig,2));
    colormap gray
    axis off
    saveas(fig, [outpath, ['/Horizontal_trig', num2str(trig), '_maxLev_angular.png']])
end
close all


% for radial modulator
for trig = 1:1:stimP.nTrig
    % for vertical
    figure('Position', [0, 0, 960, 540]);
    fig = imagesc(phaseAvgBandsLev_radial{maxLev}(:,:,trig,1));
    colormap gray
    axis off
    saveas(fig, [outpath, ['/Vertical_trig', num2str(trig), '_maxLev_radial.png']])

    figure('Position', [0, 0, 960, 540]);
    % for horizontal
    fig = imagesc(phaseAvgBandsLev_radial{maxLev}(:,:,trig,2));
    colormap gray
    axis off
    saveas(fig, [outpath, ['/Horizontal_trig', num2str(trig), '_maxLev_radial.png']])
end
close all


%% plot the Vertical minus Horizontal per modulator at the maximum level
% for angular modulator
figure('Position', [0, 0, 960, 540]);
fig = imagesc(oriDiff_angular{maxLev});
colormap gray
axis off
saveas(fig, [outpath, '/VminsH_maxLev_angular.png'])

% for radial modulator
figure('Position', [0, 0, 960, 540]);
fig = imagesc(oriDiff_radial{maxLev});
colormap gray
axis off
saveas(fig, [outpath, '/VminsH_maxLev_radial.png'])



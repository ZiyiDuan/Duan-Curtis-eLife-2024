% Create the final stimuli screen for model inputs
% The final stimuli screen includes: the target stimulus, fixation, 
% grey mask, and the full screen, the same configuration as the experiment.
% author: Zoe Duan
% date: 08/08/2023

clear; clc

addpath(genpath('/Applications/Psychtoolbox/'));

% define paths
mainpath = input('Input the path for the parent folder', 's');
% mainpath = '/System/Volumes/Data/d/DATC/datc/formatwm_ZD/public/';
datapath = [mainpath, 'models/'];
outpath = [mainpath, 'modelResults/'];
if ~exist(outpath, 'dir')
   mkdir(outpath)
end

% ask for system
exptSystem = 'scanRoom';
% load the screenP parameters based on the system
load([datapath, 'screenP', '_', exptSystem, '.mat']); 

% initialize parameters
exptP.system = exptSystem;
[exptP] = initParams(exptP, screenP);
% define interested orientations
oris = [0, 90, 15, 75, 135];

%% make the grey mask
[~, maskGrey] = makeApertMask(screenP, exptP);

%% make fixation, double check with drawOvalFix.m
% define image size of the fixation
image_size = floor(exptP.fixOutLineSize*screenP.pixels_per_deg_height);
% initialize a mask that cover the fixation with grey color
maskFix = ones(image_size, image_size)*0.5;
% Create a grid of coordinates
[x, y] = meshgrid(1:image_size, 1:image_size);
x = x-mean(x(:));
y = y-mean(y(:));
[th,r] = cart2pol(x,y);   % change the cartesian coor into polar coor

% define radiaus of all circles of the fixation
r_outCircle_1 = floor((exptP.fixOutLineSize/2)*screenP.pixels_per_deg_height);    % fixation outer circle1, black
r_outCircle_2 = floor(0.8*r_outCircle_1);                                        % fixation outer circle2, grey
r_innerFix = floor((exptP.fixSize/2)*screenP.pixels_per_deg_height);              % inner fixation, black

% draw circles from outer to inner
maskFix(r < r_outCircle_1) = 0;            % fixation outer circle1, black
maskFix(r < r_outCircle_2) = 0.5;          % fixation outer circle2, grey
maskFix(r < r_innerFix) = 0;               % inner fixation, black

%% load stimuli whose shape is
% (x, y, nTrig, nPhase, nOri)
load(['exptStim_angular_', exptSystem, '.mat']);
im_angular = finalStim;
load(['exptStim_radial_', exptSystem, '.mat']);
im_radial = finalStim;
clear finalStim stimP

%% create the final stimuli screen
% initialize the final stimuli screen
final_stimScreen_angular = zeros(size(maskGrey, 1), size(maskGrey, 2), exptP.nTrig, exptP.nPhase, numel(oris));
final_stimScreen_radial = final_stimScreen_angular;

% add the maskGrey, stimuli, and the fixation in sequence
for trig = 1:exptP.nTrig
    for phase = 1:exptP.nPhase
        for i = 1:numel(oris)
            ori = oris(i);
            ori_inds = ori+1;

            %% for the angular modulator
            % initialize the stimScreen with the maskGrey
            stimScreen_angular_curr = maskGrey;
            % move the stimuli to the center of the mask grey and replce these pixels
            stimScreen_angular_curr(size(maskGrey, 1)/2-size(im_angular, 1)/2+1 : size(maskGrey, 1)/2+size(im_angular, 1)/2, ...
                size(maskGrey, 2)/2-size(im_angular, 2)/2+1 : size(maskGrey, 2)/2+size(im_angular, 2)/2) = im_angular(:,:,trig,phase,ori_inds);
            % move the fixation to the center of the mask grey and replce these pixels
            stimScreen_angular_curr(size(maskGrey, 1)/2-size(maskFix, 1)/2+1 : size(maskGrey, 1)/2+size(maskFix, 1)/2, ...
                size(maskGrey, 2)/2-size(maskFix, 2)/2+1 : size(maskGrey, 2)/2+size(maskFix, 2)/2) = maskFix;
            % save the final stimuli screen
            final_stimScreen_angular(:, :, trig, phase, i) = stimScreen_angular_curr;

            %% for the radial modulator
            % initialize the stimScreen with the maskGrey
            stimScreen_radial_curr = maskGrey;
            % move the stimuli to the center of the mask grey and replce these pixels
            stimScreen_radial_curr(size(maskGrey, 1)/2-size(im_radial, 1)/2+1 : size(maskGrey, 1)/2+size(im_radial, 1)/2, ...
                size(maskGrey, 2)/2-size(im_radial, 2)/2+1 : size(maskGrey, 2)/2+size(im_radial, 2)/2) = im_radial(:,:,trig,phase,ori_inds);
            % move the fixation to the center of the 
            stimScreen_radial_curr(size(maskGrey, 1)/2-size(maskFix, 1)/2+1 : size(maskGrey, 1)/2+size(maskFix, 1)/2, ...
                size(maskGrey, 2)/2-size(maskFix, 2)/2+1 : size(maskGrey, 2)/2+size(maskFix, 2)/2) = maskFix;
            % save the final stimuli screen
            final_stimScreen_radial(:, :, trig, phase, i) = stimScreen_radial_curr;
            
        end
    end
end

% save results
fname = [outpath, 'final_stimScreen_angular_', exptSystem, '.mat'];
save(fname, 'final_stimScreen_angular', 'oris', '-v7.3');
fname = [outpath, 'final_stimScreen_radial_', exptSystem, '.mat'];
save(fname, 'final_stimScreen_radial', 'oris', '-v7.3');

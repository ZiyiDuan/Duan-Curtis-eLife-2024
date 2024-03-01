% generate model output for dots representations
% You need to download the following matlab tools and addpaths
% mrTools: https://github.com/justingardner/mrTools
% mgl: https://github.com/justingardner/mgl
% matlabPyrTools: https://github.com/LabForComputationalVision/matlabPyrTools
% stimulusVignetting: https://github.com/elimerriam/stimulusVignetting 
% Author: Zoe Duan
% Date: 02/21/2024

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

% load the screenP parameters based on the system
load([mainpath, 'screenP', '_', exptSystem, '.mat']); 

% initialize parameters
exptP.system = exptSystem;
[exptP] = initParams(exptP, screenP);

%% create dots image, whose length is equal to the stimulus's diameter 
% define parameters
oris = [15, 75, 135];
len = exptP.stim.imgSize;   % dva
rdot = 10; % pixels, the radius for the dot
% define sizes
nrow = round(len*screenP.pixels_per_deg_height)+2*rdot; % extend the boundary to include the whole dot
ncol = round(len*screenP.pixels_per_deg_height)+2*rdot;
cx=nrow/2+.5;
cy=ncol/2+.5;
radius = round((nrow-2*rdot)/2); % radius for the stimulus size where the dot should be

%% make a grey mask
[~, maskGrey] = makeApertMask(screenP, exptP);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make one dot with different orientations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for oii = 1:numel(oris)
    % initialize the image
    img{oii} = repmat(0.5,nrow, ncol);
    [x,y] = meshgrid(1:nrow, 1:ncol); % in pixels
    x = x-mean(x(:));
    y = y-mean(y(:));
    [th,r] = cart2pol(x,y);
    Coords=[x(:),y(:)];

    % find pixels around the edges of the current orientation
    theta = 90-oris(oii);
    % find the center of the dot
    xcent = radius * cosd(theta);
    ycent = -radius * sind(theta);
    % select pixels around the dot center
    pixSelect = sqrt((Coords(:,1)-xcent).^2+(Coords(:,2)-ycent).^2) <= rdot;
    pix_bool{oii} = reshape(pixSelect,nrow,ncol);
    % change logical 1 into darker
    img{oii}(pix_bool{oii}) = 0.2;
end

% check for img
% figure;
% imshow(img{1})

% add images together
% initialize variables to save image
dots_final = zeros(size(maskGrey, 1), size(maskGrey, 2), numel(oris));
for oii = 1:numel(oris)
    % initialize the stimScreen with the maskGrey
    img_curr = maskGrey;
    % move the dot img to the center of the mask grey and replce these pixels
    img_curr(size(maskGrey, 1)/2-size(img{oii}, 1)/2+1 : size(maskGrey, 1)/2+size(img{oii}, 1)/2, ...
        size(maskGrey, 2)/2-size(img{oii}, 2)/2+1 : size(maskGrey, 2)/2+size(img{oii}, 2)/2) = img{oii};
    % save the final stimuli screen
    dots_final(:, :, oii) = img_curr;
end

% check for img_final
% figure;
% imshow(dots_final(:,:,1))

% save results
fname = [mainpath, 'final_1dot_', exptSystem, '.mat'];
save(fname, 'dots_final', 'oris', '-v7.3');

%% generate model output for 1 dot
% define the model parameters
% define the # of orientation bands
numOrientations = 6; 
% compute the orientation bandwidth
% oriBandwidth = 360/numOrientations;   
% define the spatial frequency bandwidth
bandwidth = 0.5;        

%% build pryamid for dots
% find the dimension of image
dims=size(dots_final, [1,2]);
% compute the # of spatial frequency bands based on the bandwidth and image size
numLevels = maxLevel(dims,bandwidth);
% construct quad frequency filters
[freqRespsImag,freqRespsReal,modulDiff]= makeQuadFRs(dims,numLevels,numOrientations,bandwidth);

% initialize variables
sumBands_dots = zeros([dims,numel(oris)]);

% Cycle over oris
for oii = 1:length(oris)

    % find the image
    im = dots_final(:,:,oii);
    [pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
    % Cycle over spatial frequency bands
    for lev = 1:numLevels
        sumBandsLev_dots{lev}(:,:,oii) = zeros(dims);
        % Cycle over orientation bands
        for orientation = 1:numOrientations
            % extract frequency response
            thisBand = accessSteerBand(pyr,pind,numOrientations,lev,orientation);
            bandIm = abs(thisBand).^2;
            % sum across orientation channels
            sumBands_dots(:,:,oii) = squeeze(sumBands_dots(:,:,oii))+bandIm;
            sumBandsLev_dots{lev}(:,:,oii) = squeeze(sumBandsLev_dots{lev}(:,:,oii))+bandIm;
        end % of over orientation bands
    end % of cycling over levels
end  % of cycling over oris

% save results
save([outpath, '/modelOutput_1dot.mat'],'sumBands_dots','sumBandsLev_dots','numOrientations','bandwidth','dims','numLevels','oris','-v7.3');


%% plot model output results for 1 dot
% Cycle over oris
for i = 1:length(oris)
    % define the output path
    output_dir = [outpath, '/modelOutput_plots_1dot_', num2str(oris(i)), 'deg/'];
    if ~exist(output_dir, 'dir')
       mkdir(output_dir)
    end

    % plot averaged results
    figure('Position', [0, 0, 960, 540]);
    fig = imagesc(sumBands_dots(:,:,i));
    colormap gray
    axis off
    saveas(fig, [output_dir, 'allLev_dots.png'])
    
    % plot results for each level
    for lev=1:numLevels
        figure('Position', [0, 0, 960, 540]);
        fig = imagesc(sumBandsLev_dots{lev}(:,:,i));
        colormap gray
        axis off
        saveas(fig, [output_dir, 'lev', num2str(lev),'_dots.png'])
    end
    close all
end % of cycling over 






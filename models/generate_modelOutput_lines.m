% generate model output for lines
% You need to download the following matlab tools and addpaths
% mrTools: https://github.com/justingardner/mrTools
% mgl: https://github.com/justingardner/mgl
% matlabPyrTools: https://github.com/LabForComputationalVision/matlabPyrTools
% stimulusVignetting: https://github.com/elimerriam/stimulusVignetting 
% Author: Zoe Duan
% Date: 09/20/2023

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

%% create lines image, whose length is equal to the stimulus's diameter 
% define parameters
oris = [0, 90, 15, 75, 135];
len = exptP.stim.imgSize;   % dva
ori_filter_width = 1000;


%% make a grey mask
[~, maskGrey] = makeApertMask(screenP, exptP);


%% make fixation, double check with drawOvalFix.m
% % define image size of the fixation
% fix_size = floor(exptP.fixOutLineSize*screenP.pixels_per_deg_height);
% % initialize a mask that cover the fixation with grey color
% maskFix = ones(fix_size, fix_size)*0.5;
% % Create a grid of coordinates
% [x, y] = meshgrid(1:fix_size, 1:fix_size);
% x = x-mean(x(:));
% y = y-mean(y(:));
% [th,r] = cart2pol(x,y);   % change the cartesian coor into polar coor
% 
% % define radiaus of all circles of the fixation
% r_outCircle_1 = floor((exptP.fixOutLineSize/2)*screenP.pixels_per_deg_height);    % fixation outer circle1, black
% r_outCircle_2 = floor(0.8*r_outCircle_1);                                        % fixation outer circle2, grey
% r_innerFix = floor((exptP.fixSize/2)*screenP.pixels_per_deg_height);              % inner fixation, black
% 
% % draw circles from outer to inner
% maskFix(r < r_outCircle_1) = 0;            % fixation outer circle1, black
% maskFix(r < r_outCircle_2) = 0.5;          % fixation outer circle2, grey
% maskFix(r < r_innerFix) = 0;               % inner fixation, black


%% make a line with different orientations
nrow = round(len*screenP.pixels_per_deg_height); 
ncol = round(len*screenP.pixels_per_deg_height);
cx=nrow/2+.5;
cy=ncol/2+.5;

for oii = 1:numel(oris)
    % initialize the image
    img{oii} = repmat(0.5,nrow, ncol);
    [x,y] = meshgrid(1:nrow, 1:ncol); % in pixels
    x = x-mean(x(:));
    y = y-mean(y(:));
    [th,r] = cart2pol(x,y);
    % consider ori in both directions (upper and lower quadrants)
    % cuz projection considers direction while ori doesn't
    ori_flip = 180-oris(oii); % the orientation here is counterclockwise
    thetas = [ori_flip ori_flip-180];
    % select pixels based on two criteria 
    for tt = 1:numel(thetas)
        theta = thetas(tt);
        e=ones(1,nrow*ncol);
        [X,Y,Z]=ndgrid((1:nrow)-cx, (1:ncol)-cy, 0);
        Coords=[X(:),Y(:),Z(:)].';
        linedir=[cosd(theta), sind(theta), 0].';
        % criterion 1: coordinates should form acute angles
        AcuteAngle=(linedir.'*Coords >= 0);
        % criterion 2: projected distance squared is less than ori_filter_width
        DistLessThan=sum(cross(Coords,linedir*e).^2)<=ori_filter_width ;
        % select pixels based on these two criteria
        pix_bool_tmp{oii,tt} = AcuteAngle & DistLessThan;
        pix_bool{oii,tt} = reshape(pix_bool_tmp{oii,tt},nrow,ncol);
    end
    % select final pixels as long as they're selected in any direction
    pix{oii} = pix_bool{oii,1} | pix_bool{oii, 2};
    % Note: pix{1} is vertical - counterclockwise 
    % - pix{91} is horizontal 
    % - pix{180} is 1deg clockwise tilt away from vertical 

    % change logical 1 into darker
    img{oii}(pix{oii}) = 0.2;
    % change pix outside the circle to be grey
    img{oii}(r > nrow/2) = 0.5;
    % change pix inside the inner aperture to be grey
%     img{oii}(r < round(exptP.stim.innerSize*screenP.pixels_per_deg_height)/2) = 0.5;
    
end

% check for img
% figure;
% imshow(img{3})

%% add images together
% initialize variables to save image
lines_final = zeros(size(maskGrey, 1), size(maskGrey, 2), numel(oris));
for oii = 1:numel(oris)
    % initialize the stimScreen with the maskGrey
    img_curr = maskGrey;
    % move the line to the center of the mask grey and replce these pixels
    img_curr(size(maskGrey, 1)/2-size(img{oii}, 1)/2+1 : size(maskGrey, 1)/2+size(img{oii}, 1)/2, ...
        size(maskGrey, 2)/2-size(img{oii}, 2)/2+1 : size(maskGrey, 2)/2+size(img{oii}, 2)/2) = img{oii};
    % move the fixation to the center of the mask grey and replce these pixels
%     img_curr(size(maskGrey, 1)/2-size(maskFix, 1)/2+1 : size(maskGrey, 1)/2+size(maskFix, 1)/2, ...
%         size(maskGrey, 2)/2-size(maskFix, 2)/2+1 : size(maskGrey, 2)/2+size(maskFix, 2)/2) = maskFix;
    % save the final stimuli screen
    lines_final(:, :, oii) = img_curr;
end

% check for img_final
figure;
imshow(lines_final(:,:,1))

% save results
fname = [mainpath, 'final_line_', exptSystem, '.mat'];
save(fname, 'lines_final', 'oris', '-v7.3');






%% generate model output for lines
% define the model parameters
% define the # of orientation bands
numOrientations = 6; 
% compute the orientation bandwidth
% oriBandwidth = 360/numOrientations;   
% define the spatial frequency bandwidth
bandwidth = 0.5;        


%% build pryamid for lines
% find the dimension of image
dims=size(lines_final, [1,2]);
% compute the # of spatial frequency bands based on the bandwidth and image size
numLevels = maxLevel(dims,bandwidth);
% construct quad frequency filters
[freqRespsImag,freqRespsReal,modulDiff]= makeQuadFRs(dims,numLevels,numOrientations,bandwidth);

% initialize variables
sumBands_lines = zeros([dims,numel(oris)]);

% Cycle over oris
for oii = 1:length(oris)

    % find the image
    im = lines_final(:,:,oii);
    [pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
    % Cycle over spatial frequency bands
    for lev = 1:numLevels
        sumBandsLev_lines{lev}(:,:,oii) = zeros(dims);
        % Cycle over orientation bands
        for orientation = 1:numOrientations
            % extract frequency response
            thisBand = accessSteerBand(pyr,pind,numOrientations,lev,orientation);
            bandIm = abs(thisBand).^2;
            % sum across orientation channels
            sumBands_lines(:,:,oii) = squeeze(sumBands_lines(:,:,oii))+bandIm;
            sumBandsLev_lines{lev}(:,:,oii) = squeeze(sumBandsLev_lines{lev}(:,:,oii))+bandIm;
        end % of over orientation bands
    end % of cycling over levels

end  % of cycling over oris

% save results
save([outpath, '/modelOutput_lines.mat'],'sumBands_lines','sumBandsLev_lines','numOrientations','bandwidth','dims','numLevels','oris','-v7.3');


%% plot model output results for lines
% Cycle over oris
for i = 1:length(oris)
    % define the output path
    output_dir = [outpath, '/modelOutput_plots_lines_', num2str(oris(i)), 'deg/'];
    if ~exist(output_dir, 'dir')
       mkdir(output_dir)
    end

    % plot averaged results
    figure('Position', [0, 0, 960, 540]);
    fig = imagesc(sumBands_lines(:,:,i));
    colormap gray
    axis off
    saveas(fig, [output_dir, 'allLev_lines.png'])
    
    % plot results for each level
    for lev=1:numLevels
        figure('Position', [0, 0, 960, 540]);
        fig = imagesc(sumBandsLev_lines{lev}(:,:,i));
        colormap gray
        axis off
        saveas(fig, [output_dir, 'lev', num2str(lev),'_lines.png'])
    end
end % of cycling over 
close all



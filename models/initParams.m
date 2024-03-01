
function [exptP] = initParams(exptP, screenP)
%% exptP

exptP.nTrial = 12;

exptP.nModul = 2;
exptP.nSquare = 1;
exptP.nTrig = 2;
exptP.nPhase = 2;  
exptP.phase = linspace(0,2*pi,exptP.nPhase+1); 
    exptP.phase(end) = [];
exptP.modulator.angFreq = 5;
exptP.modulator.sinSize = 0.2;      % make sure they are the same as makeExptStim.m
exptP.nAllOri = 180;


% Aperture (grey mask)
exptP.inApert = 1.2;                             % inner aperture (degrees): diameter   
                                                 % which defines the fixation area

exptP.outApert = screenP.screenYpixels / ...     % outer aperture (degrees): diameter
    screenP.pixels_per_deg_height;               % which defines the edge of the grey mask

% fixation
exptP.fixSize = 0.4;                             % inner fixation circle (degrees): diameter
exptP.fixOutLineSize = 0.8;                      % outer circle of fixation circle: diameter
exptP.fixCol = 0;                                % inner fixation color, black
exptP.dimFixCol = 60;                            % dimmed fixation color, for ITI

% gabor stim, should be the same as makeExptStim
exptP.stim.sfStatic = 1;                         % spatial frequency for drifting gabor in test (cycles per degree)
exptP.stim.size = 11;                            % diameter for the gabor itself, degrees
exptP.stim.imgSize = exptP.stim.size+1;         % diameter including smoothed edge, degrees
exptP.stim.cont = 0.8;                           % contrast of the stimuli
exptP.stim.innerSize = exptP.inApert;   

% noise mask
exptP.mask.size = exptP.stim.imgSize;
exptP.mask.fre_u = exptP.stim.sfStatic; 
exptP.mask.fre_std = 0.3;
exptP.mask.ori_u = [];
exptP.mask.ori_std = [];
exptP.mask.cont = exptP.stim.cont;
exptP.mask.dur = 0.2; % duration each time
exptP.mask.durTimes = 3; % how many times


% orientation settings
%exptP.oriRef = [15 75 135]; 
exptP.nOriRef = 3; %6; 
exptP.oriRef = 15 + 180/exptP.nOriRef * [0:exptP.nOriRef-1];        % [15, 75, 135]
exptP.jittLim = 7;                          % +- jitterLimit degrees from the reference angle

% error & feedback
% define feedback text position
exptP.errTextCoords = [screenP.xCenter, screenP.yCenter-...
    ceil(exptP.fixSize*screenP.pixels_per_deg_height),...
    screenP.yCenter+2*ceil(exptP.fixSize*screenP.pixels_per_deg_height)];
% exptP.errorRange = 4:8:28;                     % error smaller than exptP.errorRange(1), earn most points; more than exptP.errorRange(end), earn least/no points 
% exptP.pointsRange = 0:25:25*length(exptP.errorRange);

% continuous points range, from 0-100, 2 points each step
exptP.errorRange = 0:50;
exptP.pointsRange = 100:-2:0;


% time durations (in the order of each event within a trial)
exptP.TR = 0.75;
exptP.fixDur = exptP.TR;                           
exptP.pres = exptP.TR*3 - exptP.fixDur;                               
exptP.delay = exptP.TR*16; 
exptP.respDur = exptP.TR*6;
exptP.errDur = exptP.TR*2;
exptP.iti = exptP.TR*[8 12 16];
exptP.trialDur = exptP.fixDur + exptP.pres + exptP.delay + exptP.respDur + exptP.errDur + exptP.iti;
exptP.exptDur = sum(exptP.trialDur) * exptP.nTrial/numel(exptP.trialDur); % duration of total experiment
exptP.nFramesPres = round(sec2frames(1/screenP.ifi, exptP.pres));    % number of frames for exptP.pres
exptP.nFramesResp = round(sec2frames(1/screenP.ifi, exptP.respDur)); % number of frames for exptP.respDur


return
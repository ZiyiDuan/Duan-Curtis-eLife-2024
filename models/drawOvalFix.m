
function drawOvalFix(fixSize, fixOutLineSize, fixCol, inApertSize, screenP)
    
    % define the size of the inner aperture
    discSizePix = [0 0 inApertSize*screenP.pixels_per_deg_width inApertSize*screenP.pixels_per_deg_height];
    discRectPix = CenterRectOnPoint(discSizePix, screenP.xCenter, screenP.yCenter);
    
    % define the size of the fixation outer circle1
    fixOutLineSizePix = [0 0  fixOutLineSize * screenP.pixels_per_deg_width fixOutLineSize * screenP.pixels_per_deg_height];
    fixOutLineRectPix = CenterRectOnPoint(fixOutLineSizePix, screenP.xCenter, screenP.yCenter);

    % define the size of the fixation outer circle2
    fixOutLineSizePix_2 = 0.8*[0 0  fixOutLineSize * screenP.pixels_per_deg_width fixOutLineSize * screenP.pixels_per_deg_height];
    fixOutLineRectPix_2 = CenterRectOnPoint(fixOutLineSizePix_2, screenP.xCenter, screenP.yCenter);
    
    % define the size of the inner fixation 
    fixSizePix =[0 0 fixSize * screenP.pixels_per_deg_width fixSize * screenP.pixels_per_deg_height];
    fixRectPix = CenterRectOnPoint(fixSizePix, screenP.xCenter, screenP.yCenter);

    % draw the oval fixation
    inApertCol = screenP.grey;
    Screen('FillOval', screenP.win, [inApertCol inApertCol inApertCol], discRectPix); % inner aperture, grey
    Screen('FillOval', screenP.win, [fixCol fixCol fixCol], fixOutLineRectPix); % fixation outer circle1, black
    Screen('FillOval', screenP.win, [inApertCol inApertCol inApertCol], fixOutLineRectPix_2);   % fixation outer circle2, grey
    Screen('FillOval', screenP.win, [fixCol fixCol fixCol], fixRectPix);    % inner fixation, black
    
   
return

function finalStim = makeExptStim()
    % finalStim shape for each modulator
    % (n_pix, n_pix, n_trig, n_phase, n_oris)

    % ask for system
    exptSystem = input('system: scanRoom, exptRoom, or others?', 's');

    load(['screenP', '_', exptSystem, '.mat']); 
        % -this is on apocalypse/behavioral room
        % -therefore, the parameters are not matched to the actual ones in
        % scanner room (e.g., stim size in dva)
        % -when reporting, make sure to convert the pixel size back to dva
        % in terms of the screen parameters in the scanner room
    
    %% Parameters

    % define orientations 
    oris = deg2rad([0:1:179]); 
    
    % -- these parameters are not the updated version
    % image size, should be the same as exptP, check it out!!
    stimP.size = 11;                      % diameter for the gabor itself, degrees
    stimP.imgSize = stimP.size+1;         % diameter including smoothed edge, degrees
    stimP.innerSize = 1.2;                % diameter for the oval-fixation, matchs exptP.inApert, degrees
    
    
    % conditions
    stimP.nModul = 2;                     % # of modulators; 1:angular, 2:radial
    stimP.nSquare = 1; %2;                    % # of squares; 1:square, 2:nonSquare
    stimP.nTrig = 2;                      % # of trig; 1:cos, 2:sin
    stimP.nPhase = 2;                     % # of phases within each trig
    stimP.phase = linspace(0,2*pi,stimP.nPhase+1);     % equi-distant phases
        stimP.phase(end) = [];
%     stimP.nAllOri = 180;


    % stimulus parameters
    stimP.sfStatic = 1;                   % spatial frequency for carrier (cycles per degree)
    stimP.angFreq = 5;                    % spatial frequency for modulator
    stimP.sinSize = 0.2;                  % defines edges of the aperture of the modulator
    stimP.cont = 0.8;                       % the contrast of carrier


%     stimP.square = 0;                     % 1:square wave(edges), 0:sin


    %% Prepare 
    
    visiblesize = [stimP.imgSize stimP.imgSize]*screenP.pixels_per_deg_width; % in pixels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To avoid log(0) when creating the radial modulator %
    % make sure the visiblesize is even!!! 
    % Added by Zoe, 08/07/2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    visiblesize_1 = floor(visiblesize(1));
    visiblesize_2 = ceil(visiblesize(1));
    if rem(visiblesize_1, 2) == 0
        visiblesize = floor(visiblesize);
    elseif rem(visiblesize_2, 2) == 0
        visiblesize = ceil(visiblesize);
    end

    [x,y]=meshgrid(1:visiblesize, 1:visiblesize); % in pixels
    x = x-mean(x(:));
    y = y-mean(y(:));
    
    x = Scale(x)*stimP.imgSize - stimP.imgSize/2; % in degrees
    y = Scale(y)*stimP.imgSize - stimP.imgSize/2; % in degrees
    [th,r] = cart2pol(x,y);
    
    
    %% Make modulator
    
    % Make the aperture of the modulator, which is an annulus
    % the outer circle defines the image size, stimP.size
    % the inner circle is for fixation, stimP.innerSize
    patch = double(r <= stimP.size/2); % in degrees
    patch = patch & (r > stimP.innerSize/2); % in degrees
    
    sinFilter = sin(linspace(0,pi,stimP.sinSize*screenP.pixels_per_deg_width)).^4; % in pixels
    sinFilter = sinFilter'*sinFilter;
    sinFilter = sinFilter/sum(sinFilter(:));
    
    patch = conv2(patch, sinFilter,'same');
    

    % Make the mask of the modulator
    fr = stimP.angFreq;
    a = ((4*fr+pi)/(4*fr-pi))^(2/pi);
    
    modulator_orig(:,:,1) = fr.*th; % angular
    modulator_orig(:,:,2) = log(r)/log(a); % radial
    % check if there is any INF in the matrix
    assert(~any(isinf(modulator_orig(:))), "WARNING: there is Inf in your matrix!!!")
    

    for modul = 1:stimP.nModul
        for trig = 1:stimP.nTrig
            % nonsquare wave
            if trig == 1 %cos
                modulator_tmp = cos(modulator_orig(:,:,modul));
            elseif trig == 2 %sin
                modulator_tmp = sin(modulator_orig(:,:,modul));
            end

            % square wave
            modulator_tmp(abs(modulator_tmp)<sin(pi/4)) = 0;
            modulator_tmp = sign(modulator_tmp);


            % get the mask
            mask(:,:,modul,trig) = modulator_tmp;
        end
    end

    % Make the final modulator
    modulator = mask .* patch;
    modulator_angular = squeeze(modulator(:,:,1,:));
    modulator_radial = squeeze(modulator(:,:,2,:));
    

    %% Make final stimulus for angular modul
    % initialize variables
    carrier = nan([size(modulator_angular), stimP.nPhase, numel(oris)]);
    finalStim = nan(size(carrier));
            
    % generate 
    for trig = 1:stimP.nTrig
        for phase = 1:stimP.nPhase
            for iori = 1:numel(oris)
                ori = oris(iori);
                nx = x*cos(ori) + y*sin(ori);
                % Carrier is defined by sfStatic, phase, and cont
                carrier(:,:,trig,phase,iori) = cos(nx*stimP.sfStatic*2*pi+stimP.phase(phase)) *stimP.cont*.5;
                % finalStim = Carrier .* modulator
                finalStim(:,:,trig,phase,iori) = carrier(:,:,trig,phase,iori).*...
                    modulator_angular(:,:,trig); % range between -.5 and +.5
                % make the values range between 0 and 1
                finalStim(:,:,trig,phase,iori) = 0.5 + finalStim(:,:,trig,phase,iori);
            end
        end
    end

    
    imshow(finalStim(:,:,1,1,1));

    % save results
    fname = ['exptStim_angular_', exptSystem, '.mat'];
    save(fname, 'finalStim', 'stimP', '-v7.3');


   %% Make final stimulus for radial modul
    % initialize variables
    carrier = nan([size(modulator_radial), stimP.nPhase, numel(oris)]);
    finalStim = nan(size(carrier));
            
    % generate 
    for trig = 1:stimP.nTrig
        for phase = 1:stimP.nPhase
            for iori = 1:numel(oris)
                ori = oris(iori);
                nx = x*cos(ori) + y*sin(ori);
                % Carrier is defined by sfStatic, phase, and cont
                carrier(:,:,trig,phase,iori) = cos(nx*stimP.sfStatic*2*pi+stimP.phase(phase)) *stimP.cont*.5;
                % finalStim = Carrier .* modulator
                finalStim(:,:,trig,phase,iori) = carrier(:,:,trig,phase,iori).*...
                    modulator_radial(:,:,trig); % range between -.5 and +.5
                % make the values range between 0 and 1
                finalStim(:,:,trig,phase,iori) = 0.5 + finalStim(:,:,trig,phase,iori);
            end
        end
    end

    
    imshow(finalStim(:,:,1,1,16));

    % save results
    fname = ['exptStim_radial_', exptSystem, '.mat'];
    save(fname, 'finalStim', 'stimP', '-v7.3');

return
    
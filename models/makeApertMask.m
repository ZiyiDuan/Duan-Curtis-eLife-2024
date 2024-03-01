
function [mask, maskGrey] = makeApertMask(screenP, exptP)
        % create a grey aperture mask in the screen
        
        % the outer aperture mask spanning the visual field
        
        mask = zeros(screenP.screenYpixels, screenP.screenXpixels);
        [x,y] = meshgrid(1:size(mask,2), 1:size(mask,1)); % in pixels
        x = x-mean(x(:));
        y = y-mean(y(:));
        [th,r] = cart2pol(x,y);   % change the cartesian coor into polar coor
        mask(r < exptP.outApert/2 * screenP.pixels_per_deg_height) = 1;   % change the aperture into 1(white)

        % make mask=1(white) pixels into mask=0.5(grey)
        maskGrey = mask;
        maskGrey(mask==1)=0.5;

return
function plot_recon_sep(recon, para, clim_z, figi, roi_name, savePath, varargin)
% plot spatial reconstruction maps separately for each orientation
% condition

% INPUT:
% - recon: spatial reconstruction maps for a certain ROI across subjects 
%           (n_epochs, n_moduls, n_subjs), for a certain data, 
%           the shape is (n_xgauss, n_yguass, n_conds)
% - para: parameters for plotting
% - clim_z   : min and max values for plotting
% - figi: the index of current figure
% - roi_name: the name for the plot of current roi
% - savePath: file path to save the figure
% - fit_form:
%   'fitline_all': fit weighted lines based on all pixels in the image
%   'fitline_stim': fit weighted lines based on pixels within a square 
%                   that cover the stimulus location only

%% set defaults
if isempty(clim_z)
   clim_z = [-1.2 1.2]; % when z-scoring each subj, for beta
end

if ismember('fitline_all', varargin) || ismember('fitline_stim', varargin)
    % threshold: find top 10%
    perc = para.perc;
    disp(['threshold: ', num2str(perc), ' %'])
end

% define modulators
moduls = {'radial', 'angular'};

%% parameters
n_epochs = size(recon,1);
n_moduls = size(recon,2);
n_subjs = size(recon,3);
n_conds = size(recon{1},3);

% find the radius of the stimuli
radius = para.stimr/para.binunit;
% find the center of the image
cent = repmat(round(size([-para.ecclim:para.binunit:para.ecclim],2)/2), [1 2]);
% find the xs and ys for stimuli
[stimx, stimy] = plot_stim(radius, cent); 
              
%% plot recon
figure(figi)
fig = gcf; fig.Color = 'w'; ax = gca; ax.FontSize = 14;
set(gcf,'Position',para.mapSize)

for e = 1:n_epochs
    for m = 1:n_moduls
        for jj = 1: n_conds
            % normalize the recon across subjects by using zscore
            sub_ind = 1;
            for sub = 1:n_subjs
                if ~isempty(recon{e,m,sub}) 
                    dat = recon{e,m,sub}(:,:,jj);
                    recon_z{sub_ind}(:,:,jj) = (dat-mean(mean(dat)))/std(reshape(dat, [numel(dat), 1]));
                    % update sub_ind only when it has data
                    sub_ind = sub_ind + 1;
                end
            end
            % average across subjects for each condition
            r = [cellfun(@(x) {x(:,:,jj)}, recon_z)]';
            rcat = cat(3,r{:});
            img = mean(rcat,3);
            
            %% extract the largest circle inside the image square to plot
            imgsz = [size(img,2), size(img,1)];
            % Calculate the radius and center of the square
            img_radius = imgsz(1)/2;
            center_x = imgsz(1)/2;
            center_y = imgsz(2)/2;
            % Iterate through the matrix and keep the original values inside the circle
            for x = 1:imgsz(1)
                for y = 1:imgsz(2)
                    distance = sqrt((x - center_x)^2 + (y - center_y)^2);
                    if distance > img_radius
                        img(x, y) = NaN; % Set pixel value to NaN for outside the circle
                    end
                end
            end
            
            %% start plotting 
            % find current subplot
            map_ind = (e-1)*n_moduls*n_conds + (m-1)*n_conds + jj;
            f = subplot(n_epochs*n_moduls, n_conds, map_ind);
            hold on

            % plot spatial maps
            imagesc(img,clim_z); colormap magma; hold on;
%             colorbar; 

            % plot stimuli
            plot(stimx, stimy, 'Color', [0 0 0], 'LineWidth', para.lineWidth);

            % plot fit line
            if ismember( 'fitline_all', varargin)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find the row and column for pixels that have above
                % threshold image intensity
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                imgsz = [size(img,2), size(img,1)];
                % find the image size
                dat = reshape(img, [imgsz(1)*imgsz(2) 1]);
                % reshape 2D image into 1D column
                thres = prctile(dat, 100-perc);
                % find the threshold image intensity
                ind = find(dat>thres);
                % find the pixel inds that above threshold
                [r, c] = ind2sub(size(img), ind);
                % find the corresponding row and column of inds
                % row --> y-axis, column --> x-axis

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % using fitlm: allows both weights and constraints (to go
                % through center)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [~,w]=sort(dat(ind));
                % sort the selected pixels based on image intensity in
                % ascending order and get the sort index
                % weight the linear regression model by its rank
                tmp = fitlm(c-cent(2),r-cent(1),'linear', 'Intercept', false, 'Weights', w);
                % centralize all data points to the center,
                % no intercept can constraint it go through center

                sl = tmp.Coefficients.Estimate;
                % get the estimated slope

                plot(([0:imgsz(2)]-cent(2))*sl+cent(1), 'Color', [0 0 0], 'LineWidth', para.lineWidth); 
                % plot the regression line

                %%  calculate estimated angle
                polang = cart2pol(ones(size(sl)), -sl);
                % y-axis is top-to-down, use -sl to calculate polar angle
                degang = rad2deg(polang);
                % make sure the vertical axis is angle=0 as the experiment
                calcang = 90 - degang;
                % make sure it is btw 0-180
                if calcang < 0
                    calcang = calcang + 360;
                end
%                 f.Title.String = [epoch_name, '-', moduls{m}, newline, num2str(round(calcang,2))];
%                 f.Title.FontSize = 14;

            elseif ismember( 'fitline_stim', varargin)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % constrain the line fitting to pixels that within a square
                % that cover the stimulus locations only
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find the length difference btw the image size and the
                % stimulus-square size
                diff = round([size(img,2), size(img,1)]./2) - radius;
%                 diff = diff - 10;
                % find the inds for pixels inside the stimulus-square
                image_stimSquare = img(1+diff(1):diff(1)+2*radius, 1+diff(2):diff(2)+2*radius);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find the row and column for pixels that have above
                % threshold image intensity
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                imgsz_stimSquare = [size(image_stimSquare,2), size(image_stimSquare,1)];
                % find the size for the stimulus-Square
                cent_stimSquare = round(size(image_stimSquare)/2);
                % find the center for the stimulus-Square
                dat = reshape(image_stimSquare, [imgsz_stimSquare(1)*imgsz_stimSquare(2) 1]);
                % reshape 2D image into 1D column
                thres = prctile(dat, 100-perc);
                % find the threshold image intensity
                ind = find(dat>thres);
                % find the pixel inds that above threshold
                [r, c] = ind2sub(size(image_stimSquare), ind);
                % find the corresponding row and column of inds
                % row --> y-axis, column --> x-axis

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % using fitlm: allows both weights and constraints (to go
                % through center)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [~,w]=sort(dat(ind));
                % sort the selected pixels based on image intensity in
                % ascending order and get the sort index
                % weight the linear regression model by its rank
                tmp = fitlm(c-cent_stimSquare(2),r-cent_stimSquare(1),'linear', 'Intercept', false, 'Weights', w);
                % centralize all data points to the center,
                % no intercept can constraint it go through center

                sl = tmp.Coefficients.Estimate;
                % get the estimated slope
                
                % plot the regression line
                % after fitting the line and get the slope, 
                % you should plot the results on the whole img
                imgsz = [size(img,2), size(img,1)];
                % find the image size
                plot(([0:imgsz(2)]-cent(2))*sl+cent(1), 'Color', [0 0 0], 'LineWidth', para.lineWidth); 
                % plot the regression line

                %%  calculate estimated angle
                polang = cart2pol(ones(size(sl)), -sl);
                % y-axis is top-to-down, use -sl to calculate polar angle
                degang = rad2deg(polang);
                % make sure the vertical axis is angle=0 as the experiment
                calcang = 90 - degang;
                % make sure it is btw 0-180
                if calcang < 0
                    calcang = calcang + 360;
                end
%                 f.Title.String = [epoch_name, '-', moduls{m}, newline, num2str(round(calcang,2))];
%                 f.Title.FontSize = 14;
            end

            %% Be caution!!!
            % The subplot will flip the y-axis from
            % top-to-down into bottom-to-up automatically!!!
            % Don't do this if you're not using subplot but just plot in
            % figure, which should be fine!
            set(gca, 'View', [0, -90]);
            xlim([0 size(img,2)]); ylim([0 size(img,1)]);

            axis square; axis off;
        end % of cycling conditions
    end % of cycling moduls
end % of cycling epochs

% set title based on roi_name
% fix the issue of presenting IPS0_IPS1 and IPS2_IPS3
if strcmp(roi_name, 'IPS0_IPS1')
    roi_name = 'IPS0/1';
elseif strcmp(roi_name, 'IPS2_IPS3')
    roi_name = 'IPS2/3';
end
suptitle(roi_name, 'FontSize', 40, 'FontWeight', 'bold');
legend('off');

saveas(fig, [savePath, '.png']);
set(fig, 'Renderer', 'painters'); % make text ediable
saveas(fig, [savePath,'.eps'], 'epsc');
pause(1)

end
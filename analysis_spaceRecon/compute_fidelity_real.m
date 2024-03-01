function [filterResp, fidelity] = compute_fidelity_real(recon, ori_labels, para, varargin)

% compute filtered responses and reconstruction fidelity based on spatial_recon.m 
% INPUT:
    % recon: spatial reconstruction maps for a certain ROI  
    %           the shape is (n_xgauss, n_yguass, n_conds)
    % ori_labels: orientation labels for different conds
    % ori_filter_width: projected distance squared limit
% OUTPUT:
    % filterResp: filtered responses after aligning across three conditions
    % fidelity: reconstruction fidelity based on filtered responses



%% parameters
n_conds = size(recon,3);

%% z-score each map and average across conditions 
for i = 1: n_conds
    dat = recon(:,:,i);
    recon_z(:,:,i) = (dat-mean(mean(dat)))/std(reshape(dat, [numel(dat), 1]));
end

%% create orientation filters
if ismember( 'full', varargin)
    nrow = size(recon, 1);
    ncol = size(recon, 2);
elseif ismember( 'stim', varargin)
    radius = para.stimr/para.binunit; % compute the size of stimulus radius
    nrow = 2*radius;
    ncol = 2*radius;
end

cx=nrow/2+.5;
cy=ncol/2+.5;

% define evenly spaced orientation from 0 to 180
% we're evaluating pixel intensities that are close to eval_oris
eval_oris = 0:1:180-1;

% select pixels along eval_oris
%%%%%%%%%%%%%%%%% 
for ori = 1:numel(eval_oris)
    % consider ori in both directions (upper and lower quadrants)
    % cuz projection considers direction while ori doesn't
    thetas = [eval_oris(ori) eval_oris(ori)-180];
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
        DistLessThan=sum(cross(Coords,linedir*e).^2)<=para.ori_filter_width ;
        % select pixels based on these two criteria
        pix_bool_tmp{ori,tt} = AcuteAngle & DistLessThan;
        pix_bool{ori,tt} = reshape(pix_bool_tmp{ori,tt},nrow,ncol);
    end
    % select final pixels as long as they're selected in any direction
    pix{ori} = pix_bool{ori,1} | pix_bool{ori, 2};
    % Note: pix{1} is vertical - counterclockwise 
    % - pix{91} is horizontal 
    % - pix{180} is 1deg clockwise tilt away from vertical 
    
    if ismember( 'stim', varargin)
        % put the line filter in the center of the whole recon image
        
        nrow_full = size(recon, 1);
        ncol_full = size(recon, 2);
        pix_all = false(nrow_full, ncol_full);
        pix_stim = pix{ori};

        % calculate the difference btw the center of the filtered line and
        % the center of the recon image
        diff = round((nrow_full-nrow)/2);
        % move the filtered line to the center of the recon image
        pix_all(diff+1:diff+nrow, diff+1:diff+ncol) = pix_stim;
        pix{ori} = pix_all;
    end

end

% % check for pix 
% check_ind = 1:30:180;
% check_ind = check_ind(2:end);
% figure;
% for ii=1:numel(check_ind)
%     subplot(1,numel(check_ind),ii);
%     imshow(pix{check_ind(ii)}); 
% end
% 
% % check for pix_bool
% check_ind1 = 1:30:180;
% check_ind1 = check_ind1(2:end);
% check_ind = repmat(check_ind1, [1 2]);
% figure;
% ff=0;
% for ii=1:numel(check_ind1)
%     for jj=1:2
%         ff=ff+1;
%         subplot(1,numel(check_ind),ff);
%         imshow(pix_bool{check_ind1(ii),jj}); 
%     end
% end

%% sum of pixel intensities
% for each expt condition, sum up (average doesn't work well) 
% pixel intensities in pix{ori} for each orientation filter
for ori = 1:numel(eval_oris)
    for jj = 1:n_conds
        tmp = recon_z(:,:,jj);
        % sum up
        pix_int(ori,jj) = sum(tmp(pix{ori}));
    end
end
% z-score for each condition
% z-score each column of each matrix 
pix_int_z = zscore(pix_int);

%% rotate so that the filter responses align across conditions
cent_real = numel(eval_oris)+1-ori_labels; % cent_real has to be at the middle (cent)
cent = numel(eval_oris)/2+1;
shift = cent_real - cent;
new_ind = [];
for jj = 1:n_conds
    new_ind(:,jj) = circshift(1:numel(eval_oris), -shift(jj));
    pix_int_z_rot(:,jj) = arrayfun(@(x) pix_int_z(new_ind(x,jj),jj), 1:size(new_ind,1));
end
% average across conditions and get filtered response
filterResp = mean(pix_int_z_rot, 2);


%% compute reconstruction fidelity
my_function = cosd(abs(eval_oris-90));
fidelity = mean(my_function.*filterResp');

return
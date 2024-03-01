function [recon] = do_recon_model(model_resp,pRFparams,voxel_inds,sigmaOpt,para)
% do spatial reconstruction with real data

% INPUT:
% - model_resp: model responses created by simulate_voxelResp.m (n_conds, n_voxels)
% - pRFparams: pRF parameters for all voxels
% - voxel_inds: voxel indices for current ROI, (n_voxels, 1)
% - sigmaOpt: 
%       'lim' (uses each voxel's fitted sigma value but max limit exists)
%       'fix' (fixed to certain value)
% - para: parameters for spatial reconstruction

% OUTPUT:
% - recon: spatial reconstruction maps, (n_x, n_y, n_conds)

% create x and y for guassian reconstruction
[xgauss,ygauss] = meshgrid([-para.ecclim:para.binunit:para.ecclim],[-para.ecclim:para.binunit:para.ecclim]);
% mrVista flips y positions in that - is up, + is down


% get pRF parameters for current ROI based on voxel_inds
ve = pRFparams(:, :, :, 2);
ve = ve(voxel_inds); 
ecc = pRFparams(:, :, :, 3);
ecc = ecc(voxel_inds);
sigma = pRFparams(:, :, :, 4);
sigma = sigma(voxel_inds);
x0 = pRFparams(:, :, :, 6);
x0 = x0(voxel_inds);
y0 = pRFparams(:, :, :, 7);
y0 = y0(voxel_inds);
pRFparams_ROI = [ve, ecc, sigma, x0, y0];

%% define data for reconstruction
% voxel selection
voxind = (pRFparams_ROI(:,1)>=para.velim) & ...
    (pRFparams_ROI(:,2)<=para.ecclim) & ...
    (pRFparams_ROI(:,3)<=para.sigmalim);
% sigma 
if ismember('lim', sigmaOpt)
    sigmadat = pRFparams_ROI(voxind,3);
end
% define n_conds
n_conds = size(model_resp, 1);

% z-score across ori conditions for each voxel
% to find each voxel's orientation preference among target orientations
tmp = [arrayfun(@(vox) {[zscore(model_resp(:,vox))]}, ...
    1:size(model_resp,2))]';
for vox = 1:size(model_resp,2)
    model_resp_z(:,vox) = tmp{vox};
end
% We only do this for model outputs but not for real fMRI data
% If we don't do this, we can't see the pattern

%% stimulus reconstruction
% initialize recon
recon = cell(1); 
% index of selected vox
tmp = find(voxind); 
for i = 1:n_conds
    % each cell of r is for each voxel
    if ismember('lim', sigmaOpt)
        r = arrayfun(@(v) {model_resp_z(i,v).*...
            gauss2d(xgauss,ygauss,sigmadat(v,1),[pRFparams_ROI(tmp(v),4:5)])}, ...
            1:numel(tmp));
    elseif ismember('fix', sigmaOpt)
        r = arrayfun(@(v) {model_resp_z(i,v).*...
            gauss2d(xgauss,ygauss,para.sigmafix,[pRFparams_ROI(tmp(v),4:5)])}, ...
            1:numel(tmp));
    end
    % turn cell into a matrix, where 3rd dimension is voxel
    rcat = cat(3, r{:});
    % sum across voxel dimension
    recon{1}(:,:,i) = sum(rcat, 3);
end

end
function [recon] = do_recon_weights(weights,pRFparams,voxel_inds,sigmaOpt,para)
% do spatial reconstruction with real data

% INPUT:
% - weights: weights from classifier, (n_conds, n_voxels)
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
n_conds = size(weights, 1);

%% stimulus reconstruction
% index of selected vox
tmp = find(voxind); 
for i = 1:n_conds
    % each cell of r is for each voxel
    if ismember('lim', sigmaOpt)
        r = arrayfun(@(v) {weights(i,v).*...
            gauss2d(xgauss,ygauss,sigmadat(v,1),[pRFparams_ROI(tmp(v),4:5)])}, ...
            1:numel(tmp));
    elseif ismember('fix', sigmaOpt)
        r = arrayfun(@(v) {weights(i,v).*...
            gauss2d(xgauss,ygauss,para.sigmafix,[pRFparams_ROI(tmp(v),4:5)])}, ...
            1:numel(tmp));
    end
    % turn cell into a matrix, where 3rd dimension is voxel
    rcat = cat(3, r{:});
    % sum across voxel dimension
    recon(:,:,i) = sum(rcat, 3);
end

end
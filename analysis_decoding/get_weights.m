function [weights] = get_weights(subj,patname,regsname,selname,class_args,varargin)
% get classification weights by using all data

% INPUT:
% - subj: data structure for current subject
% - patname: the name of the pattern that contain data, (n_voxels, n_TRs)
% - regsname: the name of the regressor that contain labels, (n_conditions, n_TRs)
% - selname: the name of the selector that contain run information, (1,n_TRs)
% - class_args: contain information about your classifier in a
% data structure

% OUTPUT:
% - weights: classification weights for all conditions, (n_conds, n_voxels)


% Author: Zoe Duan
% Date: 03/01/23

% z-scoring the data
subj = zscore_runs(subj,patname,selname);
patname_z = [patname, '_z'];
% get train patterns
trainpats = get_mat(subj,'pattern',patname_z);
% Load the regressors as labels
traintargs = get_mat(subj,'regressors',regsname);


% initialize results
results.header.clock = clock;

% set defaults
defaults.rand_state_int = sum(100*results.header.clock);
defaults.perfmet_functs = {'perfmet_maxclass'};
defaults.perfmet_args = struct([]);
defaults.postproc_funct = '';
defaults.ignore_unknowns = false;
defaults.allow_single_object = false;
args = propval(varargin,defaults);

% User-specified perfmet_args must be a cell array with a struct in
% each cell. If the user feeds in a struct, put it inside a cell array
if isstruct(args.perfmet_args)
    perfmet_args = args.perfmet_args;
    args = rmfield(args,'perfmet_args');
    for p=1:length(args.perfmet_functs)
        args.perfmet_args{p} = perfmet_args;
    end
    clear perfmet_args;
end

if args.rand_state_int ~= defaults.rand_state_int
    error('Rand state int argument doesn''t work properly yet');
end
rand('state',args.rand_state_int);

% Just in case the user only has one perfmet and fed it in as a
% string rather than cell array
if ~iscell(args.perfmet_functs) && ischar(args.perfmet_functs)
    % warning('Perfmet_functs should be a cell array, not a string - fixing');
    args.perfmet_functs = {args.perfmet_functs};
end

nPerfs = length(args.perfmet_functs);

sanity_check(class_args);

% Create a function handle for the classifier training function
train_funct_hand = str2func(class_args.train_funct_name);

% Call whichever training function
scratchpad = train_funct_hand(trainpats,traintargs,class_args);

% get weights 
weights = scratchpad.net.IW{1};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = sanity_check(class_args)

if ~isstruct(class_args)
  error('Class_args should be a struct');
end

if ~isfield(class_args,'test_funct_name') || ~isfield(class_args,'train_funct_name')
  error('Need to supply training and testing function names');
end

if ~ischar(class_args.test_funct_name) || ~ischar(class_args.train_funct_name)
  error('Training or testing function names have to be strings');
end

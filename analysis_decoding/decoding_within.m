function [results] = decoding_within(subj,patname,regsname,selname,class_args,varargin)
% Within dataset decoding by using cross-validation 
% refer to cross_validation function in princeton MVPA toolbox
% Calls the classifier multiple times, training and testing on
% different subsets of the data each time

% INPUT:
% - subj: data structure for current subject
% - patname: the name of the pattern that contain data, (n_voxels, n_TRs)
% - regsname: the name of the regressor that contain labels, (n_conditions, n_TRs)
% - selname: the name of the selector that contain run information, (1,n_TRs)
% - class_args: contain information about your classifier in a
% data structure

% OUTPUT:
% - subj: data structure for current subject
% - results: data structure for decoding results

% Author: Zoe Duan
% Date: 02/16/23

% z-scoring the data
subj = zscore_runs(subj,patname,selname);
patname_z = [patname, '_z'];

% create selector indices for cross-validation
subj = create_xvalid_indices(subj,selname);
selname_group = [selname, '_xval'];
% Get the names of the selectors for each iteration
selnames = find_group(subj,'selector',selname_group);
nIterations = length(selnames);
if ~nIterations
    error('No selector group to run cross-validation over - if you want to run cross_validation.m with just a single selector that you''ve created, then you need to turn it into a group first - see http://groups.google.com/group/mvpa-toolbox/browse_thread/thread/9c7dae2757205644');
end

% get the names of patterns for each iteration
patnames = cellstr(repmat(patname_z,[nIterations 1]));

% Load the regressors
regressors = get_mat(subj,'regressors',regsname);


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
results.header.rand_state_int = args.rand_state_int;
results.header.rand_state_vec = rand('state');

% Initialize the results structure
results.header.experiment = subj.header.experiment;
results.header.subj_id    = subj.header.id;

% Just in case the user only has one perfmet and fed it in as a
% string rather than cell array
if ~iscell(args.perfmet_functs) && ischar(args.perfmet_functs)
    % warning('Perfmet_functs should be a cell array, not a string - fixing');
    args.perfmet_functs = {args.perfmet_functs};
end

nPerfs = length(args.perfmet_functs);

sanity_check(class_args);

% Initialize store_perfs - this is going to store all the
% performance scores from each of the iterations, separately
% for each performance metric
%
% nPerfs x nIterations
store_perfs = [];

disp( sprintf('Starting %i cross-validation classification iterations - %s', ...
    nIterations,class_args.train_funct_name) );

for n=1:nIterations

    fprintf('\t%i',n);
    cur_iteration = [];

    cv_args.cur_iteration = n;
    cv_args.n_iterations = nIterations;

    % Set up the current selector 
    cur_selsname = selnames{n};
    selectors = get_mat(subj,'selector',cur_selsname);

    % Set up the current pattern
    cur_patname = patnames{n};
    cur_pats = get_mat(subj,'pattern',cur_patname);

    % Extract the training and testing indices from the selector
    train_idx = find(selectors==1);
    test_idx  = find(selectors==2);
    % check if there is unknown selector labels
    unused_idx  = find(selectors==0);
    unknown_idx = selectors;
    unknown_idx([train_idx test_idx unused_idx]) = [];
    if length(unknown_idx) & ~args.ignore_unknowns
        warning( sprintf('There are unknown selector labels in %s',cur_selsname) );
    end

    if isempty(train_idx) && isempty(test_idx)
        disp('No pats and targs timepoints for this iteration - skipping');
        continue
    end

    % Create the training patterns and targets
    trainpats  = cur_pats(:,train_idx);
    traintargs = regressors( :,train_idx);
    testpats   = cur_pats(:,test_idx);
    testtargs  = regressors( :,test_idx);

    % Create a function handle for the classifier training function
    train_funct_hand = str2func(class_args.train_funct_name);

    % Call whichever training function
    scratchpad = train_funct_hand(trainpats,traintargs,class_args,cv_args);

    % Create a function handle for the classifier testing function
    test_funct_hand = str2func(class_args.test_funct_name);

    % Call whichever testing function
    [acts scratchpad] = test_funct_hand(testpats,testtargs,scratchpad);

    % If a post-processing function has been specified,
    % call it on the acts and scratchpad.
    if ~isempty(args.postproc_funct)
        postproc_funct_hand = str2func(args.postproc_funct);
        [acts scratchpad] = postproc_funct_hand(acts,scratchpad);
    end

    % this is redundant, but it's the easiest way of
    % passing the current information to the perfmet
    scratchpad.cur_iteration = n;

    % Run all the perfmet functions on the classifier outputs
    % and store the resulting perfmet structure in a cell
    for p=1:nPerfs

        % Get the name of the perfmet function
        cur_pm_name = args.perfmet_functs{p};

        % Create a function handle to it
        cur_pm_fh = str2func(cur_pm_name);

        % Run the perfmet function and get an object back
        cur_pm = cur_pm_fh(acts,testtargs,scratchpad,args.perfmet_args{p});

        % Add the function's name to the object
        cur_pm.function_name = cur_pm_name;

        % Append this perfmet object to the array of perfmet objects,
        % only using a cell array if necessary
        if nPerfs==1
            cur_iteration.perfmet = cur_pm;
        else
            cur_iteration.perfmet{p} = cur_pm;
        end

        % Store this iteration's performance. If it's a NaN, the NANMEAN
        % call below will ignore it. Updated on 080910 to store NaNs.
        cur_iteration.perf(p) = cur_pm.perf;
        % nPerfs x nIterations
        store_perfs(p,n) = cur_iteration.perf(p);

    end

    % Display the performance for this iteration
    disp( sprintf('\t%.2f',cur_iteration.perf(p)) );

    % Book-keep the bountiful insight from this iteration
    t = now;
    cur_iteration.created.datetime  = datetime(t,'ConvertFrom','datenum');
    cur_iteration.train_idx         = train_idx;
    cur_iteration.test_idx          = test_idx;
    cur_iteration.unused_idx        = unused_idx;
    cur_iteration.unknown_idx       = unknown_idx;
    cur_iteration.acts              = acts;
    cur_iteration.scratchpad        = scratchpad;
    cur_iteration.header.history    = []; % should fill this in xxx
    cur_iteration.created.function  = 'cross_validation';
    cur_iteration.created.patname   = cur_patname;
    cur_iteration.created.regsname  = regsname;
    cur_iteration.created.selname   = cur_selsname;
    cur_iteration.train_funct_name  = class_args.train_funct_name;
    cur_iteration.test_funct_name   = class_args.test_funct_name;
    cur_iteration.args              = args;
    results.iterations(n) = cur_iteration;

end % for n nIterations


disp(' ');

% Show me the money
results.total_perf = nanmean(store_perfs,2);

mainhist = sprintf( ...
    'Cross-validation using %s and %s - got total_perfs - %s', ...
    class_args.train_funct_name,class_args.test_funct_name, ...
    num2str(results.total_perf'));

results = add_results_history(results,mainhist,true);



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


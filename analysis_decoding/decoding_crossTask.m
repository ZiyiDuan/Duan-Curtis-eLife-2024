function [results] = decoding_crossTask(subj_train, patname_train, regsname_train, selname_train,...
    voxel_inds_train, subj_test, patname_test, regsname_test, selname_test, voxel_inds_test, ...
    class_args, nIterations, varargin)
% Cross-task decoding, training on the control task and testing on the main
% task

% INPUT:
% - subj_train: data structure for training
% - patname_train: the name of the pattern that contain training data, (n_voxels, n_TRs)
% - regsname_train: the name of the regressor that contain training labels, (n_conditions, n_TRs)
% - selname_train: the name of the selector that contain run information for training data, (1,n_TRs)
% - voxel_inds_train: voxel inds for the training dataset, (n_voxels, 1)
% - subj_test: data structure for testing
% - patname_test: the name of the pattern that contain testing data, (n_voxels, n_TRs)
% - regsname_test: the name of the regressor that contain testing labels, (n_conditions, n_TRs)
% - selname_test: the name of the selector that contain run information for testing data, (1,n_TRs)
% - voxel_inds_test: voxel inds for the testing dataset, (n_voxels, 1)
% - class_args: contain information about your classifier 
% - nIterations: number of iterations for classification


% OUTPUT:
% - results: data structure for decoding results
  
% Author: Zoe Duan
% Date: 05/16/23  

% set the train data
subj_train = zscore_runs(subj_train, patname_train, selname_train); % z-scoring the data
patname_train_z = [patname_train, '_z'];
trainpats = get_mat(subj_train,'pattern',patname_train_z);
% choose trainpats based on voxel inds
trainpats = trainpats(voxel_inds_train, :);
% set the train labels
traintargs = get_mat(subj_train,'regressors',regsname_train);

% set the test data
subj_test = zscore_runs(subj_test, patname_test, selname_test); % z-scoring the data
patname_test_z = [patname_test, '_z'];
testpats = get_mat(subj_test,'pattern',patname_test_z);
% choose testpats based on voxel inds
testpats = testpats(voxel_inds_test, :);
% set the test labels
testtargs = get_mat(subj_test,'regressors',regsname_test);


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
results.header.experiment = subj_train.header.experiment;
results.header.subj_id    = subj_train.header.id;

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

disp( sprintf('Starting %i classification iterations - %s', ...
    nIterations,class_args.train_funct_name) );

for n=1:nIterations

    fprintf('\t%i',n);
    cur_iteration = [];

    cv_args.cur_iteration = n;
    cv_args.n_iterations = nIterations;

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
    cur_iteration.acts              = acts;
    cur_iteration.scratchpad        = scratchpad;
    cur_iteration.header.history    = []; % should fill this in xxx
    cur_iteration.train_funct_name  = class_args.train_funct_name;
    cur_iteration.test_funct_name   = class_args.test_funct_name;
    cur_iteration.args              = args;
    results.iterations(n) = cur_iteration;

end % for n nIterations

disp(' ');

% Show me the money
results.total_perf = nanmean(store_perfs,2);

mainhist = sprintf( ...
    'Cross dataset decoding using %s and %s - got total_perfs - %s', ...
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

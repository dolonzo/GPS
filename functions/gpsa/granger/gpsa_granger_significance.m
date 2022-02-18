
function varargout = gpsa_granger_significance(varargin)
% Measure granger significance, comparing the results to the null
% hypotheses.
%
% Author: A. Conrad Nied (conrad.logos@gmail.com)
%
% Changelog:
% 2012-12-13 GPS1.7 gpsa_granger_sigtests.m Created
% 2013-07-08 GPS1.8 Created to just gather the percentile data.
% 2013-07-10 State subject is now condition brain
% 2013-09-16 Added exceptions for older datasets

%% Input

[state, operation] = gpsa_inputs(varargin);

%% Prepare a report on the type or progress of the data

if(~isempty(strfind(operation, 't')))
    report.spec_subj = 0; % Subject specific?
    report.spec_cond = 3; % Condition specific?
end

%% Execute the process

if(~isempty(strfind(operation, 'c')))
    
    study = gpsa_parameter(state, state.study);
    condition = gpsa_parameter(state, state.condition);
    state.function = 'gpsa_granger_significance';
    tbegin = tic;
    
    state.subject = condition.cortex.brain;
    
    %% Load Data
    
    % Determine filenames
    rawfilename = gps_filename(study, condition, 'granger_analysis_nullhypo');
    outputfilename = gps_filename(study, condition, 'granger_analysis_results_now');
    
    % Load Data
    rawdata = load(rawfilename);
     
    %% Get P-Values
    
    % Get vector lengths
    data.N_ROIs = rawdata.N_ROIs;
    data.N_comp = rawdata.N_comp;
    data.N_time = rawdata.N_time;
    data.N_trials = rawdata.N_trials;
    
    % Get general data
    data.name = 'results';
    data.condition = condition;
    if(isfield(rawdata, 'sample_times'));
        data.sample_times = rawdata.sample_times;
    end
    if(isfield(rawdata, 'all_rois')); data.rois = rawdata.all_rois;
    elseif(isfield(rawdata, 'rois')); data.rois = rawdata.rois;
    end
    data.act = rawdata.data;
    data.results = rawdata.granger_results;
    data.srcs = rawdata.src_ROIs; %How is this used
    if(isfield(rawdata, 'sink_ROIs')); rawdata.snk_ROIs = rawdata.sink_ROIs; end
    data.snks = rawdata.snk_ROIs; %How is this used

	if iscell(data.srcs)
		data.srcs = union(data.srcs{1}', data.srcs{2});
		data.snks = union(data.snks{1}, data.snks{2}');
	end
    
    % Granger parameters
    data.model_order = rawdata.model_order;
    if(isfield(rawdata, 'pred_adapt')); data.pred_adapt = rawdata.pred_adapt;
    else data.pred_adapt = rawdata.W_gain;
    end
    data.file_rawgranger = rawdata.inputfilename;
    data.file_nullhypotheses = rawfilename;
    data.file_results = outputfilename;
    
    % Get the p_values for each point
    data.p_values = zeros(size(data.results));
%     p_values = squeeze(mean(...
%         repmat(data.results(data.snks, data.srcs, :), [1 1 1 data.N_comp])...
%         >= rawdata.total_control_granger, 4)); % >=NaN returns 0
    all_greater = zeros(size(data.results));
    for j = 1:data.N_comp
        iter_greater = data.results(data.snks, data.srcs, :) >= squeeze(rawdata.total_control_granger(:,:,:, j));
        all_greater(data.snks, data.srcs, :) = all_greater(data.snks, data.srcs, :) + iter_greater;
    end
    p_values = squeeze(all_greater./data.N_comp);
    data.p_values(data.snks, data.srcs, :) = p_values(data.snks, data.srcs, :); %do we need to change to NaN? How is this used in future?
    
    % Get the line of each percentile
    data.perc10 = zeros(size(data.results));
    perc10 = quantile(rawdata.total_control_granger, 0.9, 4); %quantile removes NaNs, so that's good
    data.perc10(data.snks, data.srcs, :) = perc10;
    data.perc05 = zeros(size(data.results));
    perc05 = quantile(rawdata.total_control_granger, 0.95, 4);
    data.perc05(data.snks, data.srcs, :) = perc05;
    data.perc01 = zeros(size(data.results));
    perc01 = quantile(rawdata.total_control_granger, 0.99, 4);
    data.perc01(data.snks, data.srcs, :) = perc01;
    data.perc005 = zeros(size(data.results));
    perc005 = quantile(rawdata.total_control_granger, 0.995, 4);
    data.perc005(data.snks, data.srcs, :) = perc005;
    data.perc001 = zeros(size(data.results));
    perc001 = quantile(rawdata.total_control_granger, 0.999, 4);
    data.perc001(data.snks, data.srcs, :) = perc001;
    
    %% Wrap up and save
    
    % Save
    save(data.file_results, '-struct', 'data');
    
    % Record the process
    gpsa_log(state, toc(tbegin));
    
end % If we should do the function

%% Add to the report concerning the progress

if(~isempty(strfind(operation, 'p')))
    study = gpsa_parameter(state.study);
    condition = gpsa_parameter(state.condition);
    
    % Predecessor: gpsa_granger_nullhypo
    report.ready = ~~exist(gps_filename(study, condition, 'granger_analysis_nullhypo'), 'file');
    report.progress = ~~exist(gps_filename(study, condition, 'granger_analysis_results'), 'file');
    
    report.finished = report.progress == 1;
end

%% Prepare the report and output

if(nargout == 1 && exist('report', 'var'));
    varargout{1} = report;
end

end % function

%COMMENTED OUT FUNCTION 4/28/21 ON
% function varargout = gpsa_granger_significance(varargin)
% % Measure granger significance, comparing the results to the null
% % hypotheses.
% %
% % Author: A. Conrad Nied (conrad.logos@gmail.com)
% %
% % Changelog:
% % 2012-12-13 GPS1.7 gpsa_granger_sigtests.m Created
% % 2013-07-08 GPS1.8 Created to just gather the percentile data.
% % 2013-07-10 State subject is now condition brain
% % 2013-09-16 Added exceptions for older datasets
% 
% %% Input
% 
% [state, operation] = gpsa_inputs(varargin);
% 
% %% Prepare a report on the type or progress of the data
% 
% if(~isempty(strfind(operation, 't')))
%     report.spec_subj = 0; % Subject specific?
%     report.spec_cond = 3; % Condition specific?
% end
% 
% %% Execute the process
% 
% if(~isempty(strfind(operation, 'c')))
%     
%     study = gpsa_parameter(state, state.study);
%     condition = gpsa_parameter(state, state.condition);
%     state.function = 'gpsa_granger_significance';
%     tbegin = tic;
%     
%     state.subject = condition.cortex.brain;
%     
%     %% Load Data
%     
%     % Determine filenames
%     rawfilename = gps_filename(study, condition, 'granger_analysis_nullhypo');
%     outputfilename = gps_filename(study, condition, 'granger_analysis_results_now');
%     
%     % Load Data
%     rawdata = load(rawfilename);
%     
%     %% Get P-Values
%     
%     % Get vector lengths
%     data.N_ROIs = rawdata.N_ROIs;
%     data.N_comp = rawdata.N_comp;
%     data.N_time = rawdata.N_time;
%     data.N_trials = rawdata.N_trials;
%     
%     % Get general data
%     data.name = 'results';
%     data.condition = condition;
%     if(isfield(rawdata, 'sample_times'));
%         data.sample_times = rawdata.sample_times;
%     end
%     if(isfield(rawdata, 'all_rois')); data.rois = rawdata.all_rois;
%     elseif(isfield(rawdata, 'rois')); data.rois = rawdata.rois;
%     end
%     data.act = rawdata.data;
%     data.results = rawdata.granger_results;
%     data.srcs = rawdata.src_ROIs;
%     if(isfield(rawdata, 'sink_ROIs')); rawdata.snk_ROIs = rawdata.sink_ROIs; end
%     data.snks = rawdata.snk_ROIs;
%     
%     % Granger parameters
%     data.model_order = rawdata.model_order;
%     if(isfield(rawdata, 'pred_adapt')); data.pred_adapt = rawdata.pred_adapt;
%     else data.pred_adapt = rawdata.W_gain;
%     end
%     data.file_rawgranger = rawdata.inputfilename;
%     data.file_nullhypotheses = rawfilename;
%     data.file_results = outputfilename;
%     
%     % Get the p_values for each point
%     data.p_values = zeros(size(data.results));
%     p_values = squeeze(mean(...
%         repmat(data.results(data.snks, data.srcs, :), [1 1 1 data.N_comp])...
%         >= rawdata.total_control_granger, 4));
%     data.p_values(data.snks, data.srcs, :) = p_values;
%     
%     % Get the line of each percentile
%     data.perc10 = zeros(size(data.results));
%     perc10 = quantile(rawdata.total_control_granger, 0.9, 4);
%     data.perc10(data.snks, data.srcs, :) = perc10;
%     data.perc05 = zeros(size(data.results));
%     perc05 = quantile(rawdata.total_control_granger, 0.95, 4);
%     data.perc05(data.snks, data.srcs, :) = perc05;
%     data.perc01 = zeros(size(data.results));
%     perc01 = quantile(rawdata.total_control_granger, 0.99, 4);
%     data.perc01(data.snks, data.srcs, :) = perc01;
%     data.perc005 = zeros(size(data.results));
%     perc005 = quantile(rawdata.total_control_granger, 0.995, 4);
%     data.perc005(data.snks, data.srcs, :) = perc005;
%     data.perc001 = zeros(size(data.results));
%     perc001 = quantile(rawdata.total_control_granger, 0.999, 4);
%     data.perc001(data.snks, data.srcs, :) = perc001;
%     
%     %% Wrap up and save
%     
%     % Save
%     save(data.file_results, '-struct', 'data');
%     
%     % Record the process
%     gpsa_log(state, toc(tbegin));
%     
% end % If we should do the function
% 
% %% Add to the report concerning the progress
% 
% if(~isempty(strfind(operation, 'p')))
%     study = gpsa_parameter(state.study);
%     condition = gpsa_parameter(state.condition);
%     
%     % Predecessor: gpsa_granger_nullhypo
%     report.ready = ~~exist(gps_filename(study, condition, 'granger_analysis_nullhypo'), 'file');
%     report.progress = ~~exist(gps_filename(study, condition, 'granger_analysis_results'), 'file');
%     
%     report.finished = report.progress == 1;
% end
% 
% %% Prepare the report and output
% 
% if(nargout == 1 && exist('report', 'var'));
%     varargout{1} = report;
% end
% 
% end % function
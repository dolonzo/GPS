function varargout = gpsa_granger_compute(varargin)
% Compute granger results
%
% Author: A. Conrad Nied (conrad.logos@gmail.com)
%
% Changelog:
% 2012.02.28 - Originally created as GPS1.6(-)/granger_compute.m
% 2012.05.17 - Last modified in GPS1.6(-)
% 2012.10.11 - Updated to GPS1.7 format
% 2012.10.17 - Fixed archiving
% 2012.12.05 - Imports name from gpsa_granger_filename
% 2013.01.14 - Added type clause for subject specific
% 2013.04.11 - GPS 1.8, Updated the status check to the new system
% 2013.04.25 - GPS1.8 Changed subset design to condition hierarchy
% 2013.06.28 - Works with updated granger analysis design

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
    state.function = 'gpsa_granger_compute';
    state.subject = condition.cortex.brain;
    tbegin = tic;

    % Format Data Filename
    if(~strcmp(condition.name, condition.cortex.roiset))
        roiset = sprintf('_ROIs-%s', condition.cortex.roiset);
    else
        roiset = '';
    end
    
    description = sprintf('%s%s%s', condition.name, roiset);
    
    % Load the data
    datafile.filename = gps_filename(study, condition, 'granger_analysis_input');
    datafile = load(datafile.filename);
    
    
    % Configure environment
    model_order = study.granger.model_order;
    pred_adapt = study.granger.W_gain;
    
    % Process Based on the Stream
    if(isfield(state, 'grangerstream') && strcmp(state.grangerstream, 'rs'))
        stream = '_Stream-rs';
        
        % Initial Kalman
        [sspace, W_all, residual] = rs_kalman(datafile.data, model_order, size(datafile.data, 1), size(datafile.data, 2), W_gain);
        
        % Granger
        granger_results = rs_granger(datafile.data, model_order, W_gain);
        
    elseif(isfield(state, 'grangerstream') && strcmp(state.grangerstream, 'seth'))
        stream = '_Stream-seth';
        
        % Initial Kalman
        [sspace, W_all, residual] = gps_sethgranger_kalman(datafile.data, model_order, pred_adapt);
        
        % Granger
        granger_results = gps_sethgranger_granger(datafile.data, model_order, pred_adapt);
        
    % This is the decoding analysis, with the activation series and
    % sub-ROIs. Toggled based on the decodingROIs structure not being
    % empty.
    elseif (isfield(datafile, 'decodingROIs') && ~isempty(datafile.decodingROIs))
        stream = '';
        
        % Initial Kalman
        [sspace, W_all, residual] = gps_kalman(datafile.data, model_order, pred_adapt);
        
        % Granger
        granger_results = gps_granger_decode(datafile.data, model_order, pred_adapt, datafile.decodingROIs, datafile.rois);

    else % Regular
        stream = '';
        
        % Initial Kalman
        [sspace, W_all, residual] = gps_kalman(datafile.data, model_order, pred_adapt);
        
        % Granger
        granger_results = gps_granger(datafile.data, model_order, pred_adapt);
    end
    
    % Check the folder exists
    folder = gps_filename(study, condition, 'granger_analysis_rawoutput_dir');
    if(~exist(folder, 'dir')); mkdir(folder); end
    
    % Prepare saving structure
    description = sprintf('%s%s', description, stream);
    datafile.description = description;
    if isstruct(granger_results)
        datafile.granger_results = granger_results(1).indices;
        datafile.granger_result_package = granger_results;
    else
        datafile.granger_results = granger_results;
    end
    datafile.sspace = sspace;
    datafile.W_all = W_all;
    datafile.residual = residual;
    datafile.model_order = model_order;
    datafile.pred_adapt = pred_adapt;
    datafile.date = datestr(now, 'yyyymmdd-hhMMss');
    
    datafile.filename = gps_filename(study, condition, 'granger_analysis_rawoutput_now');
    
    % Save new file
    save(datafile.filename, '-struct', 'datafile');
    
    % Record the process
    gpsa_log(state, toc(tbegin));
    
end % If we should do the function

%% Add to the report concerning the progress

if(~isempty(strfind(operation, 'p')))
    study = gpsa_parameter(state.study);
    condition = gpsa_parameter(state.condition);
    % Predecessor: granger_consolidate
    report.ready = ~~exist(gps_filename(study, condition, 'granger_analysis_input'), 'file');
    report.progress = ~~exist(gps_filename(study, condition, 'granger_analysis_rawoutput'), 'file');
    report.finished = report.progress == 1;
end

%% Prepare the report and output

if(nargout == 1 && exist('report', 'var')); 
    varargout{1} = report;
end

end % function

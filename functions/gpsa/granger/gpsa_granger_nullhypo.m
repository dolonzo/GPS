
function varargout = gpsa_granger_nullhypo(varargin)
% Generates null hypotheses to compare against the Granger results
% Modeled after Sato et al., 2009
%
% Author: A. Conrad Nied (conrad.logos@gmail.com)
%
% Changelog:
% 2011-04-30 Originally created by Ricky Sachdeva
% 2012-07-03 Last modified as GPS1.6/granger_significance
% 2012-11-30 Ported to GPS1.7
% 2012-12-05 Imports name from gpsa_granger_filename
% 2012-12-13 Uses parallel processing and ROI limitations
% 2013-01-14 Added single subject handling
% 2013-04-11 GPS 1.8, Updated the status check to the new system
% 2013-04-25 Changed subset design to condition hierarchy
% 2013-07-08 Updated to use gps_filename s
% 2013-07-09 Keeps sample times
% 2013-07-10 Subject name as condition brain

%% Input

[state, operation] = gpsa_inputs(varargin);

%% Prepare a report on the type or progress of the data

if(~isempty(strfind(operation, 't')))
    report.spec_subj = 0; % Subject specific?
    report.spec_cond = 1; % Condition specific?
end

%% Execute the process

if(~isempty(strfind(operation, 'c')))
    
    study = gpsa_parameter(state, state.study);
    condition = gpsa_parameter(state, state.condition);
    state.subject = condition.cortex.brain;
    state.function = 'gpsa_granger_nullhypo';
    tbegin = tic;
    
    % Determine filenames
    inputfilename = gps_filename(study, condition, 'granger_analysis_rawoutput');
    outputfilename = gps_filename(study, condition, 'granger_analysis_nullhypo_now');
    
    folder = gps_filename(study, condition, 'granger_analysis_nullhypo_dir');
    if(~exist(folder, 'dir')); mkdir(folder); end
    
    % Get inputs
    grangerfile = load(inputfilename);
    
    if isfield(grangerfile, 'decodingROIs') && ~isempty(grangerfile.decodingROIs)
        sspace_all = {grangerfile.granger_result_package(:).sspace};
        sspace = sspace_all{1}; %May need for parfor loop?
        residuals_all = {grangerfile.granger_result_package(:).residuals};
        residuals = residuals_all{1}; %May need for parfor loop?
        modifiedData = {grangerfile.granger_result_package(:).modifiedData};
        roinames = {grangerfile.rois(:).name};
        decodingROInames = {grangerfile.decodingROIs(:).name};
        [~, decodingROIlist, ~] = intersect(roinames, decodingROInames);
        decoders = 1;
        granger_results = grangerfile.granger_result_package(end).indices; %For purposes of saving
    else
        sspace = grangerfile.sspace;
        residuals = grangerfile.residual;
        decoders = 0;
        granger_results = grangerfile.granger_results; %#ok<NASGU> For purposes of saving
    end
    data = grangerfile.data;
    rois = grangerfile.rois;
    pred_adapt = grangerfile.pred_adapt;
    model_order = grangerfile.model_order;
    sample_times = grangerfile.sample_times; %#ok<NASGU>
    
    clear grangerfile;
    
    % Specific Parameters
    [N_trials, N_ROIs, N_time] = size(data);
    %N_comp = 2000;
    N_comp = study.granger.N_comp;
    
    % Determine the areas to focus on
    all_rois = {rois.name};
    src_selection = condition.granger.focus;
    inclusive_focus = 0;
    % inclusive_focus generates null hypotheses for all source-sink pairs
    % which involve at least one focus ROI. Set to zero to only generate
    % nulls between focus ROIs
    
    if(~isempty(src_selection))
        focus_ROIs = zeros(N_ROIs, 1);
        for i_ROI = 1:N_ROIs
            for isrc_selection = 1:length(src_selection)
                if(~isempty(strfind(all_rois{i_ROI}, src_selection{isrc_selection})))
                    focus_ROIs(i_ROI)=1;
                end
            end
        end
        if inclusive_focus
            src_ROIs{1} = find(focus_ROIs); %first set with focus ROIs as source
            src_ROIs{2} = setdiff(1:N_ROIs, src_ROIs{1}); %second set with all ROIs (except focus) as source
            snk_ROIs{1} = 1:N_ROIs;	%first set with all ROIs as sink
            snk_ROIs{2} = src_ROIs{1}; %second set with focus ROIs as sink
            N_src_set1 = length(src_ROIs{1});
            N_snk_set1 = length(snk_ROIs{1});
            N_src_set2 = length(src_ROIs{2});
            N_snk_set2 = length(snk_ROIs{2});
        else
            src_ROIs{1} = find(focus_ROIs);
            snk_ROIs{1} = find(focus_ROIs);
            src_ROIs{2} = [];
            snk_ROIs{2} = [];    
            N_src_set1 = length(src_ROIs{1});
            N_snk_set1 = length(snk_ROIs{1});
            N_src_set2 = 1; % zero length breaks parfor
            N_snk_set2 = 1; % zero length breaks parfor
        end
        N_src = max(length(src_ROIs{1}), length(src_ROIs{2}));
        N_snk = max(length(snk_ROIs{1}), length(snk_ROIs{2}));
    else
        src_ROIs{1} = 1:N_ROIs;
        snk_ROIs{1} = 1:N_ROIs;
        src_ROIs{2} = [];
        snk_ROIs{2} = [];
        N_src_set1 = length(src_ROIs{1});
        N_snk_set1 = length(snk_ROIs{1});
        N_src_set2 = 1; % zero length breaks parfor
        N_snk_set2 = 1; % zero length breaks parfor
        N_src = N_ROIs;
        N_snk = N_ROIs;
    end
    

    
    % Preallocate the huge null hypothesis matrix
%     total_control_granger = nan(N_snk, N_src, N_time, N_comp, 'single'); %change to nan to hopefully prevent problems when it is used since it will be a sparse matrix
    total_control_granger_set1 = nan(N_snk_set1, N_src_set1, N_time, N_comp, 'single');
    total_control_granger_set2 = nan(N_snk_set2, N_src_set2, N_time, N_comp, 'single');
    
    % Parallel Processing
    %N_parallel_processes = 8;
    N_parallel_processes = 16; %clive has 16 cores to processs 10/26/2020 ON
    %matlabpool(num2str(N_parallel_processes))
    parpool(N_parallel_processes);
    
    % Make a progress marking directory
    prog_dir = sprintf('%s/GPS/tmp_%s', gps_presets('studyparameters'), datestr(now, 'yyyymmdd_hhMMss'));
    mkdir(prog_dir);

    try
        
        % Compute ALOT of times
        parfor i_comp = 1:N_comp
        %for i_comp = 1:N_comp % for loop for debugging
            %     for i_comp = 1:N_comp
            tcomp = tic;
            
            i_null_set = 1;
            % Compare all ROIs...
            for i_src = 1:N_src_set1
                i_ROI_src = src_ROIs{i_null_set}(i_src); %#ok<*PFBNS>
                            % Step 5: For each subject, resample the residuals (random sampling
                % with replacement)
                residuals_resample = zeros(N_trials, N_ROIs, N_time);
%                 resid = zeros(N_trials, N_ROIs, N_time);
%                 ss = zeros(N_ROIs*modelorder, N_ROIs, N_time);
%                 dat = zeros(N_trials, N_ROIs, N_time);
%                 gci = zeros(N_ROIs, N_ROIs, N_time);
                
                if decoders
                    if any(decodingROIlist==i_ROI_src)
                        i_granger = find(decodingROIlist==i_ROI_src);
                    else
                        i_granger = length(decodingROIlist)+1;
                    end
                    %fprintf('Granger Set Number %d\n', i_granger);
                    resid = residuals_all{i_granger};
                    ss = sspace_all{i_granger};
                    dat = modifiedData{i_granger};
                    %disp(N_subROIs);
                else
                    resid = residuals;
                    ss = sspace;
                    dat = data;
                end
                
                N_subROIs = size(dat, 2) - N_ROIs + 1; 
                % Reorder within each trial?
                for i_trial = 1:size(dat, 1)
                    for i_ROI = 1:size(dat, 2)
                        N_reorder = N_time;
                        reorder = randperm(N_reorder);
                    
                        residuals_resample(i_trial, i_ROI, :) = resid(i_trial, i_ROI, reorder);
                    end % For all time
                end % For all trials
                
                % ... with all other ROIs
                % run time
                for i_snk = 1:N_snk_set1
                    j_ROI_snk = snk_ROIs{i_null_set}(i_snk);
                    
                    bstrap_time = zeros(N_trials, N_ROIs+N_subROIs-1, N_time);
                    
                    % Indices in the source space that record ROI i influece.
                    model_order_end = j_ROI_snk*model_order;
                    start_ind = model_order_end - (model_order-1);
                    
                    % Define the source space excluding ROI i (xi)
                    sspace_xi = ss;
                    sspace_xi(start_ind:model_order_end,i_ROI_src:i_ROI_src+N_subROIs-1,:) = 0;
                    
                    % Build a new time series
                    for n = (model_order+1):N_time % Length in terms of time
                        
                        % Build the time window for the current timepoint
                        H = dat(:, :, n - (model_order:-1:1));
                        H = reshape(H, N_trials, (N_ROIs+N_subROIs-1) * model_order);
                        
                        sspace_xi_n = squeeze(sspace_xi(:,:,n));
                        
                        % Generate a new dataset for this relationship
                        bstrap_time(:,:,n) = (H * sspace_xi_n);
                    end
                    
                    bstrap_time = bstrap_time + residuals_resample;
                    
                    bstrap_time(:,:,1:model_order) = dat(:,:,1:model_order);
                    
                    % Run granger on this
                    
                    control_granger = gps_granger(bstrap_time, model_order, pred_adapt, i_ROI_src, j_ROI_snk); % Specifing Source and Sink nodes
                    total_control_granger_set1(i_snk, i_src, :, i_comp) = single(control_granger(j_ROI_snk, i_ROI_src, :)); %#ok<PFOUS>
                    
                end
            end % For all ROIs

			if ~isempty(src_ROIs{2})
            i_null_set = 2;
            % Compare all ROIs...
            for i_src = 1:N_src_set2
                i_ROI_src = src_ROIs{i_null_set}(i_src); %#ok<*PFBNS>
                            % Step 5: For each subject, resample the residuals (random sampling
                % with replacement)
                residuals_resample = zeros(N_trials, N_ROIs, N_time);
%                 resid = zeros(N_trials, N_ROIs, N_time);
%                 ss = zeros(N_ROIs*modelorder, N_ROIs, N_time);
%                 dat = zeros(N_trials, N_ROIs, N_time);
%                 gci = zeros(N_ROIs, N_ROIs, N_time);
                
                if decoders
                    if any(decodingROIlist==i_ROI_src)
                        i_granger = find(decodingROIlist==i_ROI_src);
                    else
                        i_granger = length(decodingROIlist)+1;
                    end
                    %fprintf('Granger Set Number %d\n', i_granger);
                    resid = residuals_all{i_granger};
                    ss = sspace_all{i_granger};
                    dat = modifiedData{i_granger};
                else
                    resid = residuals;
                    ss = sspace;
                    dat = data;
                end
                
                N_subROIs = size(dat, 2) - N_ROIs + 1;
            
                % Reorder within each trial?
                for i_trial = 1:size(dat, 1)
                    for i_ROI = 1:size(dat, 2)
                        N_reorder = N_time;
                        reorder = randperm(N_reorder);
                    
                        residuals_resample(i_trial, i_ROI, :) = resid(i_trial, i_ROI, reorder);
                    end % For all time
                end % For all trials
                
                % ... with all other ROIs
                % run time
                for i_snk = 1:N_snk_set2
                    j_ROI_snk = snk_ROIs{i_null_set}(i_snk);
                    
                    bstrap_time = zeros(N_trials, N_ROIs, N_time);
                    
                    % Indices in the source space that record ROI i influece.
                    model_order_end = j_ROI_snk*model_order;
                    start_ind = model_order_end - (model_order-1);
                    
                    % Define the source space excluding ROI i (xi)
                    sspace_xi = ss;    
                    sspace_xi(start_ind:model_order_end,i_ROI_src:i_ROI_src+N_subROIs-1,:) = 0;
                    
                    % Build a new time series
                    for n = (model_order+1):N_time % Length in terms of time
                        
                        % Build the time window for the current timepoint
                        H = dat(:, :, n - (model_order:-1:1));
                        H = reshape(H, N_trials, (N_ROIs+N_subROIs-1) * model_order);
                        
                        sspace_xi_n = squeeze(sspace_xi(:,:,n));
                        
                        % Generate a new dataset for this relationship
                        bstrap_time(:,:,n) = (H * sspace_xi_n);
                    end
                    
                    bstrap_time = bstrap_time + residuals_resample;
                    
                    bstrap_time(:,:,1:model_order) = dat(:,:,1:model_order);
                    
                    % Run granger on this
                    
                    control_granger = gps_granger(bstrap_time, model_order, pred_adapt, i_ROI_src, j_ROI_snk); % Specifing Source and Sink nodes
                    total_control_granger_set2(i_snk, i_src, :, i_comp) = single(control_granger(j_ROI_snk, i_ROI_src, :)); %#ok<PFOUS>
                    
                end
            end % For all ROIs
            end % if there is a second set
            
            
            % Mark how long it has taken
            filename = sprintf('%s/%d.txt', prog_dir, i_comp);
            fid = fopen(filename, 'w');
            fprintf(fid, '%d', i_comp);
            fclose(fid);
            j_comp = length(dir(prog_dir)) - 2;
            
            history_comment = sprintf('Computed Grangers for Iteration% 5d of% 5d.',...
                j_comp, N_comp);
            duration = toc(tcomp);
            etaleft = duration * (N_comp - j_comp) / 3600 / N_parallel_processes;
            fprintf('%s in% 4.0f seconds, est. %.2f h remaining\n', history_comment, duration, etaleft);
            gpsa_log(state, duration, history_comment);
        end
        
        total_control_granger = nan(N_snk, N_src, N_time, N_comp, 'single');

        if ~isempty(src_ROIs{2})
            total_control_granger(snk_ROIs{1}, src_ROIs{1}, :, :) = total_control_granger_set1;
            total_control_granger(snk_ROIs{2}, src_ROIs{2}, :, :) = total_control_granger_set2;
        else
            total_control_granger = total_control_granger_set1;
        end
        clear total_control_granger_set1
        clear total_control_granger_set2
        %total_control_granger(snk_ROIs{1}, src_ROIs{1}, :, :) = total_control_granger_set1;
        %clear total_control_granger_set1
        %total_control_granger(snk_ROIs{2}, src_ROIs{2}, :, :) = total_control_granger_set2;
        %clear total_control_granger_set2
        
    catch errormsg
        % Clean up progress checker and parallel processing pool
        rmdir(prog_dir, 's');
        %matlabpool close
        delete(gcp('nocreate'))
        
        rethrow(errormsg);
    end
    
    % Clean up progress checker and parallel processing pool
    rmdir(prog_dir, 's');
    %matlabpool close
    delete(gcp('nocreate'))
    
    save(outputfilename, '-v7.3');
    
    message = sprintf('Finished Significance for %s %s %s', study.name, condition.name, state.subset);
    gps_email_user(message, message);
    
    % Record the process
    gpsa_log(state, toc(tbegin));
    
end % If we should do the function

%% Add to the report concerning the progress

if(~isempty(strfind(operation, 'p')))
    study = gpsa_parameter(state.study);
    condition = gpsa_parameter(state.condition);
    
    % Predecessor: gpsa_granger_compute
    report.ready = ~~exist(gps_filename(study, condition, 'granger_analysis_rawoutput'), 'file');
    report.progress = ~~exist(gps_filename(study, condition, 'granger_analysis_nullhypo'), 'file');
    
    report.finished = report.progress == 1;
end

%% Prepare the report and output

if(nargout == 1 && exist('report', 'var'));
    varargout{1} = report;
end

end % function
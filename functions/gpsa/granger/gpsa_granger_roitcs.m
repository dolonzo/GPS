function varargout = gpsa_granger_roitcs(varargin)
% Gets average waves for a condition
%
% Author: A. Conrad Nied (conrad.logos@gmail.com)
%
% Changelog:
% 2012-09-12 Originally created as GPS1.6/wave_corticalaves_set1rs5lS.m
% 2012-10-10 Updated to GPS1.7 format
% 2012-11-11 Writes to roiset now
% 2013-01-11 Forks to a separate file if the granger is single subject
% 2013-04-11 GPS 1.8, Updated the status check to the new system
% 2013-04-25 Changed subset design to condition hierarchy
% 2013-06-27 Revamped - doesn't make intermediate STC files anymore

%% Input

[state, operation] = gpsa_inputs(varargin);

%% Prepare a report on the type or progress of the data

if(~isempty(strfind(operation, 't')))
    report.spec_subj = 1; % Subject specific?
    report.spec_cond = 2; % Condition specific?
end

%% Execute the process

if(~isempty(strfind(operation, 'c')))
    
    study = gpsa_parameter(state.study);
    subject = gpsa_parameter(state.subject);
    condition = gpsa_parameter(state.condition);
    state.function = 'gpsa_granger_rtc';
    tbegin = tic;
    
    %% Set some parameters
    
    if(~isfield(state, 'progressfid'))
        state.progressfid = 1; % stdout
    end
    
    % Set this flag to 1 to see an empty figure.
    flag_image = 0; % tsg turned off 1/22
    
    if(flag_image)
        imdir = gps_filename(study, subject, 'mne_images_dir');
        if(~exist(imdir, 'dir')); mkdir(imdir); end
    end
    
    %% Load Inverse operator
    
    fprintf(state.progressfid, '\tGetting Inverse Operator\n');
    
    % Get the inverse operator
    [~, inv] = evalc('mne_read_inverse_operator(subject.mne.invfile)'); %#ok<NASGU>
    
    % Prepares other information for the inverse operator
    % Uses dSPM
    snr = 5;
    lambda2 = 1 / (snr * snr); %#ok<NASGU>
    [~, inv] = evalc('mne_prepare_inverse_operator(inv, inv.nave, lambda2, 1)');
    
    %% Get sensor data
    
    fprintf(state.progressfid, '\tGetting Sensor Data\n');
    
    events = condition.event.code;
    
    % For each event that is part of the condition
    for i_event = 1:length(events)
        event = events(i_event);
        
        % Load the event file
        eventfilename = gps_filename(study, subject, 'meg_evoked_event', ['event=' num2str(event)]);
        if(exist(eventfilename, 'file'))
            eventfile = load(eventfilename);
            
            % Add its data to the matrix of all sensor data
            if(~exist('sensordata', 'var'))
                sensordata = cat(3, eventfile.event_data.epoch);
            else
                sensordata = cat(3, sensordata, eventfile.event_data.epoch);
            end
        end % if the file exists
        
    end % for each event
    
    % If we didn't find any event data, exit, making and empty file
    % placeholder
    if(~exist('eventfile', 'var'))
        
        filename = gps_filename(study, condition, subject, 'granger_waves_rois_subject_mat');
        
        data = []; %#ok<NASGU>
        save(filename, 'study', 'subject', 'condition', 'data');
        
        % Record the process
        gpsa_log(state, toc(tbegin));
        return
    end
    
    % Select only the channels we are using
    [~, i_channels] = intersect(eventfile.channel_names, inv.eigen_fields.col_names, 'stable');
    channel_names = eventfile.channel_names(i_channels); %#ok<NASGU>
    
    % Get some data from the eventfile (technically the last one)
    sample_times = eventfile.sample_times;
    sensordata = sensordata(i_channels, :, :);
    [N_channels, N_samples, N_trials] = size(sensordata);
    
    %% Make average waves
    
    fprintf(state.progressfid, '\tMaking average waves\n');
    
    % Make sensor averages for each group of trials
    switch condition.granger.trialsperwave
        case inf
            N_waves = 1;
            wavetrials = 1:N_trials;
        case 0
            N_waves = N_trials;
            wavetrials = 1:N_trials;
        otherwise
            N_waves = floor(N_trials / condition.granger.trialsperwave);
            s = RandStream('mt19937ar', 'Seed', 0); % Fix the seed
            wavetrials = randperm(s, N_trials);
            wavetrials = reshape(wavetrials(1:(N_waves * condition.granger.trialsperwave)), N_waves, []);
            fprintf('\t\t%d Trials into %d waves of %d each\n', N_trials, N_waves, condition.granger.trialsperwave);
    end
    
    sensordata_aves = zeros(N_channels, N_samples, N_waves);
    for i_wave = 1:N_waves
        % Average the trials
        sensordata_ave = squeeze(mean(sensordata(:, :, wavetrials(i_wave, :)), 3));
        
        % Baseline the data
        sample = @(time) find(sample_times >= time, 1, 'first');
        baseline_period = sample(condition.event.basestart) : sample(condition.event.basestop);
        sensordata_ave = sensordata_ave - repmat(mean(sensordata_ave(:, baseline_period), 2), 1, N_samples);
        
        sensordata_aves(:, :, i_wave) = sensordata_ave;
    end % for each wave
    
    % Clean up the large variable
    clear sensordata
    
    %% Map to the Cortex
    
    fprintf(state.progressfid, '\tMap onto Cortex\n');
    
    % Map data onto the cortex
    N_sources = inv.nsource;
    cortexdata_waves = zeros(N_sources, N_samples, N_waves);
    for i_wave = 1:N_waves
        cortexdata_waves(:, :, i_wave) = gps_sensor2cortex(squeeze(sensordata_aves(:, :, i_wave)), inv);
    end
    
    %% Get ROI data
    
    % Get the time presets
    timestart = condition.event.focusstart/1000;
    timestop = condition.event.focusstop/1000;
    
    filename = gps_filename(state, study, condition, 'granger_rois_set_subject_mat');
    
    fprintf(state.progressfid, '\tGet ROI data\n');
    rois = load(filename);
    
    % This if presents two alternatives.  One is the 'avatar' 
    % approach, where the time series that represents each roi
    % is the most active vertex within that ROI, while the 
    % other uses a time series made up of averages of the
    % vertices within the ROI.           tsg 1/22
    if (1)
        % Use the average approach. The rois.averageActivation
        % was calculated in gpsa_granger_rois.
        N_rois = size(rois.rois, 2);
        for i_roi = 1:length(rois.rois)
            roidata(i_roi, :) = mean(cortexdata_waves([rois.rois(i_roi).vertexMap], :, :), 1);
        end
    else
        % Use the avatar approach.  This is the 'original' GPS
        % method.
        roidata = cortexdata_waves([rois.rois.decIndex], :, :);
    end

    % Also check to see if there are decoding results.  If there are,
    % we want to create time series averages for each of the exploded
    % pieces of the ROI, as well as a time series of the average error
    % for each ROI.
    decodeROIs = [];
    for iRoi = 1:length(rois.rois)
      roiDir = sprintf("%s/%s", ...
                       gps_filename(study, subject, condition, ...
                                    'decoding_analysis_subject_roi_labels_dir'), ...
                       rois.rois(iRoi).name);
      labelFiles = dir(roiDir);
      if (~isempty(labelFiles))
        % There are results for this ROI, collect all the labels for
        % each subdivision.

        % First, ignore the '.' and '..' returned by dir().
        labelFiles = labelFiles(3:end);

        % The ROI child structure names the ROI and contains a list of the 
        % sub-ROIs that it contains.
        decodeROI.name = rois.rois(iRoi).name;
        decodeROI.subrois = [];

        % Loop through the subdivisions and grab data when there is a match.
        for iLabelFile = 1:length(labelFiles)
          labelFileContents = mne_read_label_file(sprintf("%s/%s", roiDir, ...
                                                          labelFiles(iLabelFile).name));
          subroi.parent = rois.rois(iRoi).name;
          % Grab the name of this sub-roi from the file name.
          [~, subroi.name, ~] = fileparts(labelFiles(iLabelFile).name);
          subroi.vertexMap = [];
          subroi.vertices = [];

          %% Check to make sure at least one of the subROI's vertices
          %% has an activation time series associated with it.
          if (sum(ismember(labelFileContents.vertices, rois.rois(iRoi).vertices)) > 0)
            % If so, find it.
            for iVertex = 1:length(labelFileContents.vertices)
              if (ismember(labelFileContents.vertices(iVertex), ...
                           rois.rois(iRoi).vertices))
                 i = find(labelFileContents.vertices(iVertex) == ...
                          rois.rois(iRoi).vertices);
                 subroi.vertices = [subroi.vertices rois.rois(iRoi).vertices(i)];
                 subroi.vertexMap = [subroi.vertexMap rois.rois(iRoi).vertexMap(i)];
              end
            end
            % Grab the activation data for each of the vertices in this subROI 
            % that have any.
            subroi.activationData = mean(cortexdata_waves([subroi.vertexMap], :), 1);
          else
            subroi.activationData = [];
            disp(sprintf("Missing a time series for %s", ...
                         labelFiles(iLabelFile).name));
          end
          decodeROI.subrois = [decodeROI.subrois subroi];
        end
        % While we're here, also grab the error time series...
         
        resultFile = sprintf("%s/%s/%s_%s*.mat", ...
          gps_filename(study, subject, condition, 'decoding_analysis_subject_results_dir'),...
          rois.rois(iRoi).name, subject.name, rois.rois(iRoi).name);
        resultFile = dir(resultFile);
        avgAccuracy = load(fullfile(resultFile.folder, resultFile.name));
        % ... and pack it into the decodingROI struct.
        names = fieldnames(avgAccuracy);
        ind = contains(names, 'accuracy');
        decodeROI.avgAccuracy = avgAccuracy.(names{ind});
          
        % DS: just store the accuracy time series in addition to the
        % accuracy values. Consolidate interpolates anyway, so if we
        % provide the accuracy time series, we don't need to worry about
        % sampling issues/mismatched time windows (as long as the zero
        % point remains the same.
        decodeROI.accTimes = avgAccuracy.Time;
        
        decodeROIs = [decodeROIs, decodeROI];
      end
    end
        
    if(flag_image)
        % Format subject name
        subjname = subject.name;
        underscores = strfind(subjname, '_');
        for uloc = fliplr(underscores)
            subjname = [subjname(1:uloc - 1) '\' subjname(uloc:end)];
        end
        % Format condition name
        subsname = condition.name;
        underscores = strfind(subsname, '_');
        for uloc = fliplr(underscores)
            subsname = [subsname(1:uloc - 1) '\' subsname(uloc:end)];
        end
        
        figure(4)
        clf
        hold on;
        selected_rois = [];
        selected_roinames = {};
        for i_ROI = 1:length(rois.rois)
            name = rois.rois(i_ROI).name;
            if(sum(strcmp({'STG', 'MTG', 'SMG'}, rois.rois(i_ROI).area)))
                selected_rois(end + 1) = i_ROI; %#ok<AGROW>
                underscores = strfind(name, '_');
                for uloc = fliplr(underscores)
                    name = [name(1:uloc - 1) '\' name(uloc:end)];
                end
                selected_roinames = cat(1, selected_roinames, name);
            end
        end
        %plot(sample_times, squeeze(mean(roidata(selected_rois, :, :), 3)));
        %xlim([timestart-.1, timestop]);
        %xlabel('Time (s)');
        %ylabel('Activation (Am)');
        %titlestr = sprintf('%s %s Selected ROI Waves', subjname, subsname);
        %title(titlestr);
        %legend(selected_roinames)
        
        filename = sprintf('%s/%s_%s_act_selectrois.png', imdir, subject.name, condition.name);
        saveas(gcf, filename);
    end % If we are making images
    
    % Save ROI Data
    folder = gps_filename(study, condition, 'granger_waves_rois_dir');
    if(~exist(folder, 'dir')); mkdir(folder); end
    
    filename = gps_filename(study, condition, subject, 'granger_waves_rois_subject_mat');
    
    savedata.data = roidata;
    savedata.rois = rois.rois;
    savedata.decodeROIs = decodeROIs;
    savedata.sample_times = sample_times;
    savedata.study = study.name;
    savedata.subject = subject.name;
    savedata.condition = condition.name; %#ok<STRNU>
    save(filename, '-struct', 'savedata');
    
    % Record the process
    gpsa_log(state, toc(tbegin));
    
end % If we should do the function

%% Add to the report concerning the progress

if(~isempty(strfind(operation, 'p')))
    study = gpsa_parameter(state.study);
    subject = gpsa_parameter(state.subject);
    condition = gpsa_parameter(state.condition);
    
    if(sum(strcmp(condition.subjects, subject.name)))
        report.ready = mean([~~exist(gps_filename(state, study, condition, 'granger_rois_set_subject_mat'), 'file'), ...
            length(dir(gps_filename(study, subject, 'meg_evoked_gen'))) >= 1]);
        report.progress = ~~exist(gps_filename(study, condition, subject, 'granger_waves_rois_subject_mat'), 'file');
    else
        report.ready = 0; report.progress = 0; report.applicable = 0;
    end
    report.finished = report.progress == 1;
end

%% Prepare the report and output

if(nargout == 1 && exist('report', 'var'));
    varargout{1} = report;
end

end % function


    

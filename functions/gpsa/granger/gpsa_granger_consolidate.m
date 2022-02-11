function varargout = gpsa_granger_consolidate(varargin)
% Assembles Granger data for the condition
%
% Author: A. Conrad Nied
%
% Changelog:
% 2012.02.24 - Originally created as GPS1.6(-)/granger_assemble.m
% 2012.05.17 - Last modified in GPS1.6(-)
% 2012.10.11 - Updated to GPS1.7 format
% 2012.11.11 - stoff
% 2012.12.05 - Imports name from gpsa_granger_filename
% 2013.01.11 - Added clause for single subject analysis
% 2013.01.14 - Finished single subject additions
% 2013.04.25 - GPS1.8 Changed subset design to condition hierarchy
% 2013.06.28 - Works with updated granger analysis design
% 2013.07.10 - Fixed single subject processing, gets subject list from
% condition

%% Input

[state, operation] = gpsa_inputs(varargin);

%% Prepare a report on the type or progress of the data

if(~isempty(strfind(operation, 't')))
    report.spec_subj = 0; % Subject specific?
    report.spec_cond = 3; % Condition specific?
end

%% Execute the process

if(~isempty(strfind(operation, 'c')))
     
    study = gpsa_parameter(state.study);
    condition = gpsa_parameter(state.condition);
    state.function = 'gpsa_granger_consolidate';
    tbegin = tic;
    
    subjects = condition.subjects;
    subjects_final = subjects;
    N_subjects = length(subjects);
    
    % Get the time presets
    timestart = condition.event.focusstart;
    timestop = condition.event.focusstop;
    
    % For all subjects
    for i_subject = 1:N_subjects
        subject = gpsa_parameter(subjects{i_subject});
        fprintf('\tProcessing subject: %s', subject.name);
        
        % Load the ROI information for the subject condition
        roifilename = gps_filename(study, condition, subject, 'granger_waves_rois_subject_mat');
        roidata = load(roifilename);
        
        if(~isempty(roidata.data))
            % If this is the first subject, allocate the matrices
            if(i_subject == 1 || ~exist('sample_times', 'var'))
                start = min(0, timestart);
                %start = timestart - 100; %10/29/20 ON wanted to have a shorter time start above zero. 
                stop = min(roidata.sample_times(end)*1000, timestop + 100);
                sample_times = (start:stop) / 1000;
                N_samples = length(sample_times);
                
                rois = roidata.rois;
                N_ROIs = length(rois);
                
                N_waves = 0;
                data = zeros(N_waves, N_ROIs, N_samples);

                decodingROIs = [];
            end % if this is the first subject
            
            if(isempty(roidata.data))
                subjects = setdiff(subjects_final, i_subject);
            else
                [N_ROIs_subject, ~, N_waves_subject] = size(roidata.data);
                if(N_ROIs_subject ~= N_ROIs);
                    error('Inconsistent amount of ROIs between the subject and previous');
                end
                for i_wave = 1:N_waves_subject
                    for i_ROI = 1:N_ROIs
                        data(i_wave + N_waves, i_ROI, :) = interp1(roidata.sample_times, squeeze(roidata.data(i_ROI, :, i_wave)), sample_times);
                    end

		    % What is happening here is that we are rearranging the
		    % sub-roi data so that it is sorted by ROI and subject,
		    % instead of subject and ROI, which is how it arrives.  It
		    % arrives as a struct containing a collection of subject
		    % data, each of which contains a set of ROIs, each with
		    % accuracy data and a collection of subROI data. We leave
		    % with a collection of ROIs, with the subject data combined
		    % into big blocks.

                    % For this subject (wave) make an array of the subROI
                    % time series collection and the avgAccuracy series,
                    % too.
                    for (i_ROI = 1:length(roidata.decodeROIs))
                      decodingROIs(i_ROI).name = roidata.decodeROIs(i_ROI).name;
                      decodingROIs(i_ROI).avgAccuracy(i_subject, :) = ...
                          interp1(roidata.decodeROIs(i_ROI).accTimes, roidata.decodeROIs(i_ROI).avgAccuracy, sample_times);
                      N_subROIs = length(roidata.decodeROIs(i_ROI).subrois);
                      for (i_subROI = 1:N_subROIs)
                        decodingROIs(i_ROI).subROIdata(i_subject, i_subROI, :) = ...
                          interp1(roidata.sample_times, roidata.decodeROIs(i_ROI).subrois(i_subROI).activationData, sample_times);
                      end
                    end
                end
                fprintf(' %d waves saved.\n', N_waves_subject);
                N_waves = N_waves + N_waves_subject;
            end
        end
    end % for each subject
    
    % Make sure the input directory exists
    folder = gps_filename(study, condition, 'granger_analysis_input_dir');
    if(~exist(folder, 'dir')); mkdir(folder); end
    
    % Save the file to input
    savefile.filename = gps_filename(study, condition, 'granger_analysis_input');
    savefile.subjects = subjects_final;
    savefile.N_subjects = length(subjects_final);
    savefile.N_waves = N_waves;
    savefile.N_ROIs = N_ROIs;
    savefile.N_samples = N_samples;
    savefile.data = data;
    savefile.sample_times = sample_times;
    savefile.rois = rois;
    savefile.decodingROIs = decodingROIs;
    save(savefile.filename, '-struct', 'savefile');
    % Record the process
    gpsa_log(state, toc(tbegin));
    
end % If we should do the function

%% Add to the report concerning the progress

if(~isempty(strfind(operation, 'p')))
    study = gpsa_parameter(state.study);
    condition = gpsa_parameter(state.condition);
    
    % Predecessor: granger_roitcs
    report.ready = 0;
    
    for i_subject = 1:length(condition.subjects);
        subject = gpsa_parameter(condition.subjects{i_subject});
        
        roifilename = gps_filename(study, condition, subject, 'granger_waves_rois_subject_mat');
        report.ready = report.ready + double(~~exist(roifilename, 'file')) / length(condition.subjects);
    end % for each subject
    report.progress = ~~exist(gps_filename(study, condition, 'granger_analysis_input'), 'file');
    report.finished = report.progress == 1;
end

%% Prepare the report and output

if(nargout == 1 && exist('report', 'var'));
    varargout{1} = report;
end

end % function

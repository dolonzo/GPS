function varargout = gpsa_granger_rois(varargin)
% Template function for all other functions
%
% Author: A. Conrad Nied (conrad.logos@gmail.com)
%
% Changelog:
% 2012.01.04 - Originally created as GPS1.6(-)/granger_vertex.m
% 2012.08.27 - Also based on GPS1.6(-)/wave_gather_rois.m
% 2012.09.13 - Last modified in GPS1.6(-)
% 2012.10.09 - Updated to GPS1.7 format
% 2012.10.12 - Images the ROIs now
% 2012.10.17 - Fixed problem in max_i, it was indexed ascending not
% descending so this program found the weakest point in the ROI...
% 2012.11.11 - Takes subset.cortex.roiset instead of .roidir
% 2013.04.11 - GPS 1.8, Updated the status check to the new system
% 2013.04.25 - Changed subset design to condition hierarchy
% 2013.06.25 - Fixed some bugs in the new GPS system and added annotation
% section.
% 2013.07.10 - Performed using the condition brain rather than a study one
% 2019.01-03 - Added explicit pathname references to environment vars.  -tsg

%% Input

[state, operation] = gpsa_inputs(varargin);

%% Prepare a report on the type or progress of the data

if(~isempty(strfind(operation, 't')))
    report.spec_subj = 1; % Subject specific?
    report.spec_cond = 3; % Condition specific?
end

%% Execute the process

if(~isempty(strfind(operation, 'c')))
    
    study = gpsa_parameter(state.study);
    condition = gpsa_parameter(state.condition);
    if(strcmp(state.subject, condition.cortex.brain));
        subject.name = condition.cortex.brain;
        subject.type = 'subject';
        subject.mri.dir = sprintf('%s/%s', study.mri.dir, subject.name);
    else
        subject = gpsa_parameter(state.subject);
    end
    state.function = 'gpsa_granger_rois';
    tbegin = tic;
    
    % Configure
    % Set this to 1 to see an empty figure. Not sure why this is.
    flag_image = 0; % tsg turned off 1/22 
    
    % Import their brain
    brain = gps_brain_get(subject);
    
    flag_morph = 1;
    hemis = {'lh'; 'rh'};
    N_hemis = 2;
    
    % Get the ROI directories
    roidir = gps_filename(study, condition, 'granger_rois_set_dir');
    subject_roidir = gps_filename(state, study, condition, 'granger_rois_set_subject_dir');
    
    if(~isempty(dir([roidir '/*label'])))
        
        % Get the time presets
        timestart = condition.event.focusstart;
        timestop =  condition.event.focusstop;
        
        % If necessary, rewrite the labels for the ROI to be mapped on the
        % subject's brain
        if (flag_morph)
            % Remove old morphed maps
            labels = gps_filename(state, study, condition, 'granger_rois_set_subject_labels_gen');
            if(~isempty(dir(labels)))
                delete(labels)
            end
            
            tlocal = tic;
            %% Added explicit ref to mnehome.  -tsg
            unix_command = sprintf('%s $MNE_ROOT/bin/mne_morph_labels --from %s --to %s --labeldir %s --smooth 5',...
                state.setenv, condition.cortex.brain, subject.name, roidir);
            [~, ~] = unix(unix_command);
            gpsa_log(state, unix_command, toc(tlocal));
        end
        
        N_ROIs = 0;
        regions_cortical = zeros(brain.N, 1);
        activity_cortex = zeros(brain.decN, 1);
        activity_ROIs = [];
        
        % For Each Hemisphere
        for i_hem = 1:N_hemis
            hemi = hemis{i_hem};
            
            %% Get Cortical Activity
            
            if(strcmp(subject.name, condition.cortex.brain))
                actfile = sprintf('%s-%s.stc', gps_filename(study, condition, 'mne_stc_avesubj'), hemi);
            else
                actfile = sprintf('%s-%s.stc', gps_filename(subject, condition, 'mne_stc'), hemi);
            end
            
            if strcmp(study.name, 'MPS1')
                folder = sprintf('%s/%s/stcs', study.meg.dir, state.subject);
                segment = condition.name(3:end-1);
                actfile = sprintf('%s/%s_%s_median_mne-%s.stc', folder, state.subject, segment, hemi);
            end
            
            if(~exist(actfile, 'file'))
                [actfile, path] = uigetfile(actfile(1:find(actfile == '/', 1, 'last')-1));
                actfile = sprintf('%s/%s', path, actfile);
            end
            actdata = mne_read_stc_file1(actfile);
            
            sample_times = actdata.tmin + ((1:size(actdata.data, 2)) - 1) * actdata.tstep;
            start = max(find(sample_times > timestart/1000, 1, 'first') - 1, 1);
            stop = min(find(sample_times < timestop/1000, 1, 'last') + 1, length(sample_times));
            time_focus = start:stop;
            
            if(strcmp(hemi, 'lh'))
                activity_cortex(1:brain.decN_L) = mean(actdata.data(:, time_focus), 2);
            else
                activity_cortex(brain.decN_L + 1 : end) = mean(actdata.data(:, time_focus), 2);
            end
            
            % Find all of the ROIs
            labels = sprintf('%s/*%s.label', subject_roidir, hemi);
            labels = dir(labels);
            
            N_hemiROIs = length(labels);
            rois_files = cell(N_hemiROIs, 1);
            ROIcents = zeros(N_hemiROIs, 1);
            
            % For each label, find the point with the highest PLV activation and
            % make this the representative of the ROI
            for i_hemiROI = 1:N_hemiROIs
                roi_file = labels(i_hemiROI).name;
                rois_files{i_hemiROI} = roi_file;
                
                % Locate the vertex with the highest PLV
                labelfile = sprintf('%s/%s', subject_roidir, roi_file);
                labelfile_contents = mne_read_label_file(labelfile);
                [activity, ~, vertices] = mne_label_time_courses1(labelfile, actfile); % mid is time
                % Changed to update the label vertex list to be 1 based
                % (matching the stc and inverse operators)  dos 3/22
                
                
                % Avatar method. We characterize an ROI by using the
                % activation time series from a single vertex within
                % it. We choose the one with the "highest" activation
                % value, defined as the integral of activation over
                % the time range.                            tsg 1/22
                [~, max_i] = sort(sum(activity(:, time_focus), 2));

                % Do not repeat vertices (this is a quick not the best fix)
                while(sum(vertices(max_i(end)) == ROIcents))
                    max_i(end) = [];
                end
                centroid_vertex = vertices(max_i(end));
                ROIcents(i_hemiROI) = centroid_vertex;
                
                %% Prepare structure
                
                % Name
                roi.name = roi_file(1 : end - 9);
                roi.file = labelfile;
                
                % Features
                dashes = strfind(roi.name, '-');
                if(isempty(dashes))
                    dashes = strfind(roi.name, '_');
                end
                roi.hemi = roi.name(1);
                roi.area = roi.name(dashes(1) + 1 : dashes(2) - 1);
                roi.num = str2double(roi.name(dashes(2) + 1: end));
                
                % Vertex
                if(strcmp(hemi, 'lh'))
                    roi.centroid = centroid_vertex;
                else
                    roi.centroid = centroid_vertex + double(brain.N_L);
                end
                
                % Decimated Vertex
                roi.decIndex = find(brain.decIndices == roi.centroid);
                
                if(isempty(roi.decIndex));
                    error('Did not find the decimated index of the %s', roi.name);
                end

                % Also add a mapping between vertex numbers and the
                % corresponding activation time series, also for use
                % later.                                      tsg 1/22
                roi.vertexMap = [];
                for kVertex = 1:length(vertices)
                  if (strcmp(hemi, 'lh'))
                    vertexNumber = vertices(kVertex);
                  else
                    vertexNumber = vertices(kVertex) + double(brain.N_L);
                  end
                  roi.vertexMap(kVertex) = find(brain.decIndices == vertexNumber);
                end
                
                % FS Parcellation
                roi.aparcI = brain.aparcI(roi.centroid);
                roi.aparc = brain.aparcShort{roi.aparcI};
                
                % Medial or Lateral?
                switch roi.area
                    case {'STS', 'cMFG', 'AG', 'ITG', 'LOC', 'LOrb', 'MTG', 'Oper',...
                            'ParsOper', 'ParsTri', 'Tri', 'OrbIFG', 'Calc',...
                            'postCG', 'preCG', 'rMFG', 'SFG', 'SPC', 'STG', 'SMG',...
                            'FPol', 'TPol', 'Aud', 'Insula', 'ParsOrb'}
                        roi.side = 'L';
                    case {'ParaHip','Medial', 'caCing', 'Ent', 'CC', 'Isth', 'Cun',...
                            'ParaC', 'MOrb', 'Fusi', 'pCing', 'preCun', 'raCing', 'Ling'}
                        roi.side = 'M';
                    otherwise
                        roi.side = 'M';
                end
                if(roi.hemi == 'L')
                    roi.vertices = vertices;
                    roi.allvertices = labelfile_contents.vertices + 1;
                else
                    roi.vertices = vertices + brain.N_L;
                    roi.allvertices = labelfile_contents.vertices + 1 + brain.N_L;
                end
                
                % Add to structure
                N_ROIs = N_ROIs + 1;
                rois(N_ROIs) = roi; %#ok<AGROW>
                
                % Put in regions
                if(hemi == 'L')
                    regions_cortical(vertices + 1) = N_ROIs;
                else
                    regions_cortical(vertices + brain.N_L + 1) = N_ROIs;
                end
                
                % Activity
                activity_ROIs(N_ROIs, :) = activity(max_i(1), time_focus); %#ok<AGROW>
            end
            
            % Save old version
            subjectlabelmat = sprintf('%s/%s_ROIcentroids-%s.mat',...
                roidir, subject.name, hemi);
            labelfilessave = rois_files; %#ok<NASGU>
            save(subjectlabelmat, 'ROIcents', 'labelfilessave');
            
        end % For Each Hemisphere
        
        % Make annotation
        filename = sprintf('%s/label/lh.%s.annot', subject.mri.dir, condition.name);
        if(state.override || ~exist(filename, 'file'))
            gps_labels2annot(state, subject_roidir);
        end
        gpsa_mri_brain2mat(state);
        brain = gps_brain_get(state);
        
        % Display on the cortex
        if(flag_image)
            figure(2)
            clf
            set(gcf, 'Units', 'Pixels');
            set(gcf, 'Position', [10, 10, 800, 600]);
            set(gca, 'Units', 'Normalized');
            set(gca, 'Position', [0, 0, 1, 1]);
            
            % Set parameters
            drawdata = brain;
            drawdata.act.data = activity_cortex;
            drawdata.act.p = [80 90 95];
            drawdata.regions = regions_cortical;
            drawdata.points = rois;
            
            options.overlays.name = 'act';
            options.overlays.percentiled = 'p';
            options.overlays.decimated = 1;
            options.overlays.coloring = 'h';
            options.shading = 1;
            options.curvature = 'bin';
            options.sides = {'ll', 'rl', 'lm', 'rm'};
            options.fig = gcf;
            options.axes = gca;
            options.centroids = 1;
            options.centroids_color = ones(length(rois), 3) * 0.5;
            options.parcellation = condition.name;
            options.parcellation_text = 1;
            options.parcellation_overlay = 'top';
            options.parcellation_border = 2;
            
            % Draw mean activity
            %%%%%%%%%%%%%%% Commented out by SA 2015-06-16, because crashed
            %%%%%%%%%%%%%%% with error: java.lang.IllegalArgumentException:
            %%%%%%%%%%%%%%% adding a container to a container on a different GraphicsDevice 
            %gps_brain_draw(drawdata, options);
            
            frame = getframe(gcf);
            filename = sprintf('%s/%s_%s_rois.png',...
                subject_roidir, subject.name, condition.name);
            imwrite(frame.cdata, filename);
            
            clear drawdata options;
        end
        
        % Summarize into a structure
        roistruct.rois = rois;
        roistruct.N = length(rois);
        roistruct.N_L = sum([rois.hemi] == 'L');
        
        %% Save
        filename = gps_filename(state, study, condition, 'granger_rois_set_subject_mat');
        save(filename, '-struct', 'roistruct');
    end % If the subset has ROIs
    
    % Record the process
    gpsa_log(state, toc(tbegin));
    
    % If this is the last subject, do it again for the average subject if
    % possible
    if(strcmp(subject.name, condition.subjects{end}) && ~strcmp(subject.name, condition.cortex.brain))
        state.subject = condition.cortex.brain;
        gpsa_granger_rois(state)
    end
    
end % If we should do the function

%% Add to the report concerning the progress

if(~isempty(strfind(operation, 'p')))
    study = gpsa_parameter(state.study);
    subject = gpsa_parameter(state.subject);
    condition = gpsa_parameter(state.condition);
    
    % Predecessor: GPSr and the forward solution (for brain2mat)
    if(sum(strcmp(condition.subjects, subject.name)))
        roidir = gps_filename(study, condition, 'granger_rois_set_labels_gen');
        
        report.ready = (double(~isempty(dir(roidir))) + ...
            double(~~exist(subject.mne.fwdfile, 'file'))) / 2;
        report.progress = ~~exist(gps_filename(state, study, condition, 'granger_rois_set_subject_mat'), 'file');
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


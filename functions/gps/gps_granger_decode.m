function gciPkg = gps_granger_decode(data, model_order, pred_adapt, decodingROIs, rois, varargin)
% Performs Granger analysis on the data
%
% Adapted from the original gps_granger to do some tricky stuff with ROIs associated
% with the decoding data.
%
% Changelog:
% 2022-02-03 - Copied from the original.
%
% Input: The data matrix of sensors by time, the model order, and the
% factor that the prediction error adapts by, and more specifications if
% they are useful at limiting the data overhead
% Output: A set of computed granger causality indices, one for each decoding ROI.

%% Setup

[N_subjects, N_ROIs, N_time] = size(data);

% Get the *names* of the ROIs in the data.
roiList = {};

for i_ROI = 1:length(rois)
  roiList{end + 1} = rois(i_ROI).name;
end

% Source and Receiving ROIs to be looked at in Granger
if (nargin > 5)
  src_ROI_names = varargin{1};
  rcv_ROI_names = varargin{2};
  src_ROIs = varargin{1};
  rcv_ROIs = varargin{2};
else % Do all
  src_ROIs = 1:N_ROIs; % ROIs that are being tested for influence from
  rcv_ROIs = 1:N_ROIs; % ROIs that are being tested for influence on
end

%% Analysis
% Gather the index numbers of the decoding ROIs.
decodingROIlist = [];
for i_decodeROI = 1:length(decodingROIs)
  decodingROIlist = [decodingROIlist, find(string(decodingROIs(i_decodeROI).name) == roiList)];
end

% Save these so they can be reused each time through the loop.
all_src_ROIs = src_ROIs;
all_rcv_ROIs = rcv_ROIs;

% A loop across decoding ROIs. One additional iteration without exploding
% any ROIs but replacing all decoders with accuracy data
for i_decodeROI = 1:length(decodingROIs)+1 
  % Create the granger_causality_indices array that will receive the results for
  % this ROI.
  granger_causality_indices = zeros(N_ROIs, N_ROIs, N_time);
  % Reset src and rcv ROIs
  src_ROIs = all_src_ROIs;
  rcv_ROIs = all_rcv_ROIs;
    
  % This is the index into the ROI list of this decoding ROI.
  if i_decodeROI > length(decodingROIs) %if the last iteration with no exploding
      sourceROIindex = 0; %zero sourceROIindex
      sourceROIdata = []; %emtpy sourceROIdata
      N_subROIs = 1; %we always use N_subROIs-1, so set N_subROIs to 1 to cancel out
  else
      sourceROIindex = decodingROIlist(i_decodeROI);

      % Grab the activation time series for each of the subROIs.
      sourceROIdata = decodingROIs(i_decodeROI).subROIdata;

      % Count the subROIs.
      [~, N_subROIs, ~] = size(sourceROIdata);
  end
  
  % Substitute the accuracy time series for all the "good" ROIs that are
  % not the one we're looking at just now.
  modifiedData = data;
  for i_goodROI = 1:length(decodingROIlist)
    if (decodingROIlist(i_goodROI) == sourceROIindex)
      continue;  % Skip the one we're looking at just now.
    end

    % Define a scale factor here because the avgAccuracy data is many
    % orders of magnitude different than the activation data.  Not
    % sure that matters, but aesthetically it seems wrong.
    scale = 1.0e-10;
  
    % Substitute the (scaled) avgAccuracy time series for the activation
    % time series in the block of time series data.  The procedure is
    % different if this is the first or last ROI in the list.
    if (decodingROIlist(i_goodROI) == 1) % Are we first?
      modifiedData = cat(2, scale*reshape(decodingROIs(i_goodROI).avgAccuracy, ...
        N_subjects, 1, N_time), modifiedData(:, 2:N_ROIs, :));
    elseif (decodingROIlist(i_goodROI) == N_ROIs) % Are we last?
      modifiedData = cat(2, modifiedData(:, 1:(N_ROIs-1), :), ...
          scale*reshape(decodingROIs(i_goodROI).avgAccuracy, N_subjects, 1, N_time));
    els % Middle, apparently.
      targetROI = decodingROIlist(i_goodROI);
      modifiedData = cat(2, modifiedData(:, 1:(targetROI-1), :), ...
          scale*reshape(decodingROIs(i_goodROI).avgAccuracy, N_subjects, 1, N_time),...
          modifiedData(:, (targetROI+1):N_ROIs, :));
    end
  end

  % For the ROI we're looking at just now, stick all of its subROI data in
  % place of the single time series for sourceROIindex.  Again, this looks a 
  % little different if it's the first or last ROI in the list.
  if (sourceROIindex == 1) % Are we first?
    modifiedData = cat(2, sourceROIdata, modifiedData(:, 2:N_ROIs, :));
  elseif (sourceROIindex == N_ROIs) % Are we last?
    modifiedData = cat(2, modifiedData(:, 1:(N_ROIs-1), :), sourceROIdata);
  elseif (sourceROIindex ~= 0) % Middle, apparently and not zero.
    modifiedData = cat(2, modifiedData(:, 1:(sourceROIindex-1), :),...
        sourceROIdata, modifiedData(:, (sourceROIindex+1):N_ROIs, :));
  end
  
  % Change the indices for sink and source ROIs based on the additional
  % subROIs.
  adj_srcROIs = [src_ROIs(src_ROIs<=sourceROIindex), ... 
                 src_ROIs(src_ROIs>sourceROIindex)+N_subROIs-1];
  adj_rcvROIs = [rcv_ROIs(rcv_ROIs<=sourceROIindex), ...
                 rcv_ROIs(rcv_ROIs>sourceROIindex)+N_subROIs-1];
  
  % Exclude accuracy time series from source ROIs
  [~, nonaccIndices] = setdiff(src_ROIs, setdiff(decodingROIlist, sourceROIindex));
  src_ROIs = src_ROIs(nonaccIndices);
  adj_srcROIs = adj_srcROIs(nonaccIndices);
  
  % Exclude exploded time series from sink ROIs
  [~, nonexpIndices] = setdiff(rcv_ROIs, sourceROIindex);
  rcv_ROIs = rcv_ROIs(nonexpIndices);
  adj_rcvROIs = adj_rcvROIs(nonexpIndices);
  
  prediction_error_standard = gps_kalman(modifiedData, model_order, pred_adapt);
  % Excluded ROI x prediction_error matrix
  prediction_errors_withoutROI = zeros(N_ROIs, N_ROIs+N_subROIs-1, ...
                                       N_ROIs+N_subROIs-1, N_time); 

  % For Each ROI that influences another ROI
  loopindex = 0;
  for i_ROI = src_ROIs
      
    loopindex = loopindex+1;
    i_adj_ROI = adj_srcROIs(loopindex);
      
    % Remove this srcROI
    dataReducer = ones(N_ROIs+N_subROIs-1, 1);
    if i_ROI == sourceROIindex %if it's the exploded ROI, remove all the time series
      dataReducer(i_adj_ROI:i_adj_ROI+N_subROIs-1) = 0;
    else  %... otherwise, just zero out the one time series.
      dataReducer(i_adj_ROI) = 0;
    end
    dataReducer = find(dataReducer);
    data_reduced = modifiedData(:, dataReducer, :);
      
    % Get the output from the Kalman filter. This output is smaller since
    % there are fewer input ROIs
    kalmanOutput = gps_kalman(data_reduced, model_order, pred_adapt);
    
    %Put the output into the full size array
    prediction_error_withoutROI = zeros(N_ROIs+N_subROIs-1, N_ROIs+N_subROIs-1, N_time);
    prediction_error_withoutROI(dataReducer, dataReducer, :) = kalmanOutput;
    
    % Park that output data where it belongs in the collection of output data.    
    prediction_errors_withoutROI(i_ROI, :, :, :) = prediction_error_withoutROI;
  end % For Each ROI in src_ROIs loop

  % For Each Time Point (starting at when you can because of the model order)
  for i_time = (model_order + 1):N_time;
    src_loopindex = 0;
    rcv_loopindex = 0;
    % For Each ROI (this one is the reduced one)
    for i_ROI = src_ROIs
      src_loopindex = src_loopindex+1;
      i_adj_ROI = adj_srcROIs(src_loopindex);
        
      rcv_loopindex = 0;
      % For Each ROI (Again)
      for j_ROI = rcv_ROIs
        rcv_loopindex = rcv_loopindex+1;
        j_adj_ROI = adj_rcvROIs(rcv_loopindex);
        % Get the prediction error in both models
        prediction_error_withoutROI = prediction_errors_withoutROI(i_ROI, j_adj_ROI, j_adj_ROI, i_time);
        prediction_error = prediction_error_standard(j_adj_ROI, j_adj_ROI, i_time);
            
        % Calculate the difference in the prediction error with and without the ROI
        granger_causality_index = log(prediction_error_withoutROI / prediction_error);
            
        granger_causality_indices(j_ROI, i_ROI, i_time) = granger_causality_index;
      end % For Each ROI
    end % For Each ROI
  end % For Each Timepoint

  % Tuck new granger_causality_indices into a struct with ROI labels.
  if i_decodeROI > length(decodingROIs)
      gciPkg(i_decodeROI).name = 'No Exploded';
  else
      gciPkg(i_decodeROI).name = decodingROIs(i_decodeROI).name;
  end
  gciPkg(i_decodeROI).indices = granger_causality_indices;

end % End loop over decoding ROIs.

end % function

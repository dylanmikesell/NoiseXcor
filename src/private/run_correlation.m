%--------------------------------------------------------------------------
% Main routine to run correlations
%--------------------------------------------------------------------------
function run_correlation( iDay, nDays, DataTable, CoordTable, fileType, ...
    corrParam, corrFilter, outputDir, startDay, stationTag, channel_list )

fprintf('\nCorrelating day %d of %d\n', iDay, nDays); % print information
todaysFiles = DataTable(:,iDay); % get data files for this day

noIdx = strcmp( todaysFiles, 'N'); % find stations w/o data
todaysFiles(noIdx) = []; % remove these empty rows

if ~isempty( todaysFiles )
    
    w = waveform; % allocate an empty waveform
    for iFile = 1 : numel( todaysFiles ) % load all files into single waveform
        w(iFile) = waveform( todaysFiles{iFile}, fileType );
    end
    fprintf('Finished loading %d waveforms\n', numel( todaysFiles ) );
    
    % Apply a filter prior to correlation if requested
    if corrFilter.prefilt
        w1 = taper( w, 0.01, 'cosine' );
        f  = filterobject( 'B', [corrFilter.fmin, corrFilter.fmax], 2 );
        w  = filtfilt( f, w1 );
    end
    
    % run the correlations over each time window
    correlate_windows( w, corrParam, corrFilter, outputDir, startDay, ...
        CoordTable, stationTag, channel_list );
else
    disp('No data to correlate toady.');
end

end % run_correlation()
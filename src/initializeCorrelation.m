function initializeCorrelation( dataBaseName, corrFilter, corrParam, np )
% 
% USAGE: initializeCorrelation( dataBaseName, corrFilter, corrParam, np )
% 
% DESCRIPTION: runs the actual correlations on one or multiple cores.
% 
% INPUT:
%   dataBaseName
%   corrFilter
%   corrParam
%   np
% OUTPUT:
%   NULL
% 

% Last modified: 15 May 2017


load( dataBaseName ); % will load 'stationData' structure

%--------------------------------------------------------------------------
% Check USER inputs
%--------------------------------------------------------------------------

if corrParam.overlapPercent > 1
    error('OverlapPercent cannot be larger than 1.0');
end

%--------------------------------------------------------------------------
% Check whether or not we should overwrite the output directory
%--------------------------------------------------------------------------
outputDir = fullfile( stationData.projectDirectory, ...
    'COR', num2str(corrFilter.filterNum, '%02d' ) );

[success, message, messageID] = checkOutputDir( outputDir ); % make directory

if success
    fprintf('Created OUTPUT directory: %s\n', outputDir)
else
    error('Did not create output directory. Check permissions');
end
clear success message messageID parts

%--------------------------------------------------------------------------
% Write processing parameters to output directory
%--------------------------------------------------------------------------
writeCorParamReadme( outputDir, corrParam, corrFilter, ...
    fullfile( stationData.projectDirectory, dataBaseName ) );

%--------------------------------------------------------------------------
% run through days and correlate
%--------------------------------------------------------------------------
tic % start the timer

nDays     = numel( stationData.Date );
DataTable = stationData.DataTable; % declare outside of parfor loop
fileType  = stationData.fileType; % declare outside of parfor loop
dates     = stationData.Date; % declare outside of parfor loop 

if np > 1 % do parallel computation
    
    if isempty( gcp('nocreate') ) % Start the parallel pool if needed
        parpool( np );
    end
    
    parfor(iDay = 1:nDays, np) % loop through each day and correlate
        runCorrelation( iDay, nDays, DataTable, fileType, corrParam, ...
            corrFilter, outputDir, dates(iDay) );     
    end
    
    % turn off the parallel pool after correlations finish
    if isempty( gcp('nocreate') ) == 0
        delete( gcp('nocreate') );
    end
    
else % do serial computation
    
    for iDay = 1 : nDays % loop through each day and correlate
        
        runCorrelation( iDay, nDays, DataTable, fileType, corrParam, ...
            corrFilter, outputDir,  dates(iDay) );     

    end
    
end % end parallel or serial

toc % print elasped time to the screen

end

%--------------------------------------------------------------------------
% Main routine to run correlations
%--------------------------------------------------------------------------
function runCorrelation( iDay, nDays, DataTable, fileType, corrParam, ...
    corrFilter, outputDir, startDay )

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
    % run the correlations
    correlateWindows( w, corrParam, corrFilter, outputDir, startDay );
else
    disp('No data to correlate toady.');
end

end % runCorrelation.m
%--------------------------------------------------------------------------

% Junk below here from old versions of code


%%
%
%
% inputDir = '/hammer/DATA/Llaima/junk/'; % directory of day-long MAT files
% % inputDir = '/peteroa/RAW_MATfiles/'; % directory of day-long MAT files
%
% outputDir = '/hammer/DATA/Llaima/junk2/'; % directory for output of the correlations
%
%
%
%
% % global logFileID
% % logFileID = fopen('log.txt','at');
% % fprintf(logFileID, 'text line number 1 \n');
%
%
%
%
%
%
%
%     fprintf( '\nProcessing day %d of %d\n', kk, numel(file) );
%
%     % load a test matrix
%     X = load( fullfile( file(kk).folder, file(kk).name ) ); % load to temporary variable
%     W = X.W; % get the waveforms from temporary variable
%
%     correlateWindows( W, corrParam, corrFilter, outputDir );
%
%
%
% %     % set default that whitening has NOT occurred
% %     for ii = 1 : numel(W)
% %         W(ii) = addfield( W(ii), 'isWhite', false );
% %     end
% %
% %     %----------------------------------------------------------------------
% %     % Process waveforms
% %     %----------------------------------------------------------------------
% %     if strcmp( wMethod, 'ftn' ) % apply frequency-time normalization in one step
% %
% %         FB = [wfmin, wfmax]; % frequency band for whitening
% %         W  = waveformWhiten( W, FB, wMethod ); % whiten the waveforms
% %
% %     else % apply time normalization followed by frequency normalization
% %
% %         % apply amplitude normalization
% %         if ampNorm
% %             W = waveformNormalization( W, timeNorm, wfmin );
% %         end
% %
% %         if whiten
% %             FB = [wfmin, wfmax]; % frequency band for whitening
% %             W  = waveformWhiten( W, FB, wMethod ); % whiten the waveforms
% %         end
% %
% %     end
% %     %----------------------------------------------------------------------
% %     % Correlate waveform windows
% %     %----------------------------------------------------------------------
% %     runDayCorrelations( W, windowLengthMinutes, overlapPercent, tMaxOut, smoothMethod, Wn, K, outputDir );
% %
% end
% toc % stop timer
% %

%
%
%



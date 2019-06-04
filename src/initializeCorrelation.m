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

%--------------------------------------------------------------------------
% Check USER inputs
%--------------------------------------------------------------------------
check_inputs( corrParam, dataBaseName )

load( dataBaseName, 'stationData' ); % will load 'stationData' structure

%--------------------------------------------------------------------------
% Setup output directory
%--------------------------------------------------------------------------
outputDir = setup_COR_dir( stationData, corrFilter );

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
CoordTable = stationData.CoordTable; % declare outside of parfor loop
stationTag = stationData.stationTag; % declare outside of parfor loop
fileType  = stationData.fileType; % declare outside of parfor loop
dates     = stationData.Date; % declare outside of parfor loop 
channel_list = stationData.channel_list;

if np > 1 % do parallel computation
    
    if isempty( gcp('nocreate') ) % Start the parallel pool if needed
%         parpool( np );
        parpool('local',np);
    end
    
    parfor(iDay = 1:nDays, np) % loop through each day and correlate
        runCorrelation( iDay, nDays, DataTable, CoordTable, fileType,...
            corrParam, corrFilter, outputDir,  dates(iDay), stationTag, ...
            channel_list);     
    end
    
    % turn off the parallel pool after correlations finish
    if isempty( gcp('nocreate') ) == 0
        delete( gcp('nocreate') );
    end
    
else % do serial computation
    for iDay = 1 : nDays % loop through each day and correlate
        runCorrelation( iDay, nDays, DataTable, CoordTable, fileType,...
            corrParam, corrFilter, outputDir,  dates(iDay), stationTag, ...
            channel_list);     
    end % correlation of iDay
end % end parallel or serial

toc % print elasped time to the screen

end % initializeCorrelation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subroutines/functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% check inputs
%--------------------------------------------------------------------------
function check_inputs( corrParam, dataBaseName )

% Check the station coordinate file exists
assert( isfile( dataBaseName ) == 1,...
    'Correlation database file does not exist. Check file PATH!')

if corrParam.overlapPercent > 1
    error('OverlapPercent cannot be larger than 1.0');
end

end % check_inputs()

%--------------------------------------------------------------------------
% Setup the output directory
%--------------------------------------------------------------------------
function outputDir = setup_COR_dir( stationData, corrFilter )

cor_directory = fullfile( stationData.projectDirectory, 'COR' );

if exist(cor_directory,'dir') ~= 7 % directory does not exist
    [success,message,messageID] = mkdir(cor_directory); % make new directory
    if success
        fprintf('Created COR directory: %s\n', cor_directory);
    else
        disp(messageID);
        disp(message);
        error('Did not create COR directory. Check permissions');
    end
else % directory exists
    disp('Found COR directory in the project directory');
end

% Check whether or not we should overwrite the filter output directory
outputDir = fullfile( stationData.projectDirectory, ...
    'COR', num2str(corrFilter.filterNum, '%02d' ) );

if exist( outputDir, 'dir' ) ~= 7 % the directory does not exist
    [success,message,messageID] = mkdir(outputDir); % make new directory
    if success
        fprintf('Created filter directory: %s\n', outputDir);
    else
        disp(messageID);
        disp(message);
        error('Did not create filter directory. Check permissions');
    end    
else % the directory does exist
    
    % Determine if the user wants to overwrite the directory or not
    saveQuestion = questdlg(...
        ['Found an existing filter with number ', num2str(corrFilter.filterNum, '%02d' ),'. What do you want to do?'], ...
        'Filter folder already exists', ...
        'Overwrite','Append','Cancel','Cancel');
    
    % Handle response
    switch saveQuestion
        case 'Overwrite'
            disp(['Overwriting ' outputDir]);
            
            % Remove the directory and all its contents first
            [success, message] = rmdir(outputDir,'s'); % recursively remove old directory
            if success
                disp('Successfully removed old directory.');
            else
                disp(message);
                error('Could not remove old directory. Check permissions');
            end
            % Make the new directory
            [success,message,messageID] = mkdir(outputDir); 
            if success
                fprintf('Created filter directory: %s\n', outputDir);
            else
                disp(messageID);
                disp(message);
                error('Did not create filter directory. Check permissions');
            end
        case 'Cancel'
            disp('Canceling this job.');
            return % go back to MATLAB prompt after ending processing
        case 'Append'
            disp('Appending data. Caution -- existing data may be overwritten.')
    end
end

end % setup_COR_dir()

%--------------------------------------------------------------------------
% Main routine to run correlations
%--------------------------------------------------------------------------
function runCorrelation( iDay, nDays, DataTable, CoordTable, fileType, ...
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
        w1 = taper(w,0.01,'cosine');
        f = filterobject('B',[corrFilter.fmin, corrFilter.fmax],2);
        w = filtfilt(f,w1);
    end
    
    % run the correlations
    correlateWindows( w, corrParam, corrFilter, outputDir, startDay, ...
        CoordTable, stationTag, channel_list );
else
    disp('No data to correlate toady.');
end

end % runCorrelation()

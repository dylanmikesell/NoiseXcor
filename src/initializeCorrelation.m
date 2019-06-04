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
write_cor_param_readme( outputDir, corrParam, corrFilter, ...
    fullfile( stationData.projectDirectory, dataBaseName ) );

%--------------------------------------------------------------------------
% run through days and correlate
%--------------------------------------------------------------------------
tic % start the timer

nDays        = numel( stationData.Date );
DataTable    = stationData.DataTable; % declare outside of parfor loop
CoordTable   = stationData.CoordTable; % declare outside of parfor loop
stationTag   = stationData.stationTag; % declare outside of parfor loop
fileType     = stationData.fileType; % declare outside of parfor loop
dates        = stationData.Date; % declare outside of parfor loop 
channel_list = stationData.channel_list;

if np > 1 % do parallel computation
    
    if isempty( gcp('nocreate') ) % Start the parallel pool if needed
%         parpool( np );
        parpool('local',np);
    end
    
    parfor(iDay = 1:nDays, np) % loop through each day and correlate
        run_correlation( iDay, nDays, DataTable, CoordTable, fileType,...
            corrParam, corrFilter, outputDir,  dates(iDay), stationTag, ...
            channel_list);     
    end
    
    % turn off the parallel pool after correlations finish
    if isempty( gcp('nocreate') ) == 0
        delete( gcp('nocreate') );
    end
    
else % do serial computation
    for iDay = 1 : nDays % loop through each day and correlate
        run_correlation( iDay, nDays, DataTable, CoordTable, fileType,...
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
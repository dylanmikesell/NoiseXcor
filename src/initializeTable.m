function stationData = initializeTable( projectDirectory, dataDirectory, ... 
    dataBaseName, fileType, data_structure, startDate, endDate )
%  
% USAGE: stationData = initializeTable( projectDirectory, dataDirectory, ...
%   dataBaseName, fileType, data_structure, startDate, endDate )
% 
% DESCRIPTION: This function will create a database file that contains the
% information needed to run cross correlations of seismic noise. This
% database contains the path and file names of all seismic data files
% within a specified time frame. The code will search through all
% subfolders of 'dataDirectory' and return paths to all 'fileType' data.
% 
% A database file is written in the '.mat' format and used later for
% parellel processing of raw seismic data followed by correlation.
% 
% INPUT: 
%   projectDirectory = path to folder where output files will be written
%   dataDirectory    = path to data (not down into the NET/STA/ folders)
%   dataBaseName     = name of the datebase file we create
%   fileType         = 'sac', 'seed', 'miniseed' (only 'sac' implemented)
%   data_structure   = 'SDS', 'BUD', 'IDDS', 'PDF' (only 'BUD' implemented)
%   startDate        = first day of database ['YYYY-MM-DD HH:MM:SS.FFF']
%   endDate          = last day of database ['YYYY-MM-DD HH:MM:SS.FFF']
% OUTPUT:
%  stationData = a database structure containing file path information.
%   The database file is populated and also written to disk.
% 
% NOTE: I copied data formats from MSNoise. Currently only BUD is implemented
%
% data_structure['SDS']  = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
% data_structure['BUD']  = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
% data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
% data_structure['PDF']  = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
% 
% Written by: Dylan Mikesell (dylanmikesell@boisetate.edu)
% Last modified: 22 February 2017 
% 
% Example:
% 
% projectDirectory = '/hammer/SOFTWARE/NoiseXcor/Example';
% dataDirectory    = fullfile( projectDirectory, 'DATA');
% dataBaseName     = 'Llaima2015_db.mat'; 
% fileType         = 'sac';
% data_structure   = 'BUD'; 
% startDate        = '2015-01-01 00:00:00'; 
% endDate          = '2015-04-01 00:00:00';
% 
% stationData = initializeTable( projectDirectory, dataDirectory, ...
%   dataBaseName, fileType, data_structure, startDate, endDate );


%--------------------------------------------------------------------------
% start scanning through the data files
%--------------------------------------------------------------------------

% get list of files
fileList = dir( fullfile( dataDirectory, '**/*.*') );

% remove directories from the list
dirIdx = [fileList.isdir];
fileList(dirIdx) = [];
clear dirIdx;

% get file extension and compare to the
names = {fileList.name}; % names of files found
nFiles = numel( names ); % number of files found
ext = cell(nFiles,1); % allocate vector
for ii = 1 : nFiles
    [~, ~, tmpExt] = fileparts( names{ii} ); %file extension
    ext{ii} = lower( tmpExt(2:end) ); % remove '.' in front of extension
end

fileIdx = strcmp( ext, fileType ); % determine files that match data format
fileList = fileList( fileIdx ); % keep only those files
clear ext names nFiles tmpExt fileIdx ii

% Get station and network data
stationFolders = transpose( unique( {fileList.folder} ) );
nStations = numel( stationFolders ); % number of unique NET/STA paths
stationName = cell(nStations,1); % allocate
networkName = cell(nStations,1); % allocate
for ii = 1 : nStations
    tmp = split(stationFolders{ii},'/');
    stationName{ii} = tmp(end);
    networkName{ii} = tmp(end-1);
end
[~,netIdx] = unique( [networkName{:}] );
networkName = networkName(netIdx);
clear ii tmp netIdx

%% Process the file list and populate the station table

% Build a structure to keep track of days with data
stationData.stationList = stationFolders;
stationData.Date = datenum( startDate ) : datenum( endDate );

nDays = numel( stationData.Date ); % number of days in the table

stationData.DataTable        = cell( nStations, nDays );
stationData.DataTable(:)     = {'N'}; % 'N' stands for no data
stationData.files            = fileList;
stationData.dataDirectory    = dataDirectory;
stationData.projectDirectory = projectDirectory;
stationData.fileType         = fileType;

for iFile = 1 : numel( fileList )
    % fprintf( 'Processing: %s\n', fileList(iFile).name );
    
    switch data_structure % account for different naming conventions
        case 'SDS' % data_structure['SDS'] = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
            % not implemented yet
        case 'BUD' % data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
            fileInfo     = split( fileList(iFile).name, '.' );
            fileStation  = fileInfo{1};
            fileNetwork  = fileInfo{2};
            fileLocation = fileInfo{3};
            fileChannel  = fileInfo{4};
            fileYear     = fileInfo{5};
            fileJulDay   = fileInfo{6};     
        case 'IDDS' % data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
            % not implemented yet
        case 'PDF' % data_structure['PDF'] = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
            % not implemented yet
        otherwise
            error('No data structures besides BUD implemented yet');
    end
    
    dataDate = datenum( [fileYear '.' fileJulDay], 'yyyy.dd' );
    
    if dataDate <= datenum( endDate ) && dataDate >= datenum( startDate )
        
        dayIdx = dataDate - datenum( startDate ) + 1; % ordinal number

        fprintf( 'Adding waveform at %s.%s.%s.%s on %s\n', ...
            fileNetwork, fileStation, fileLocation, fileChannel, ...
            datestr( dayIdx-1 + datenum(startDate) ) );

        thisStation = [fileNetwork '/' fileStation];
        k = strfind( stationFolders, thisStation );
        
        % old empty string test implementation
        % stationIdx = ~cellfun(@isempty,k);
        % version below is faster according to:
        % https://www.mathworks.com/matlabcentral/answers/16383-how-do-i-check-for-empty-cells-within-a-list
        stationIdx = ~cellfun('isempty', k);
        
        stationData.DataTable(stationIdx,dayIdx) = { fullfile( fileList(iFile).folder, fileList(iFile).name ) }; % Insert file with path
    end
end

% save database file
save( fullfile( projectDirectory, dataBaseName ), 'stationData');

% write some stuff to the screen for the user
if numel(networkName) > 1
    fprintf('Found %d unique stations over %d networks.\n', nStations, numel(networkName) );
else
    fprintf('Found %d unique stations in %d network.\n', nStations, numel(networkName) );
end
fprintf( 'Found %d data files.\n', numel( fileList ) );
fprintf( 'Finished writing database file: %s\n', fullfile( projectDirectory, dataBaseName ) );
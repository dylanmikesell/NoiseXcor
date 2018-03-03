function stationData = initializeTable( project_directory, data_directory, ...
    database_name, file_type, data_structure, start_date, end_date, channel_list )
%
% USAGE: stationData = initializeTable( project_directory, data_directory, ...
%   database_name, file_type, data_structure, start_date, end_date, channel_list )
%
% DESCRIPTION: This function will create a database file that contains the
% information needed to run cross correlations of seismic noise. This
% database contains the path and file names of all seismic data files
% within a specified time frame. The code will search through all
% subfolders of 'data_directory' and return paths to all 'file_type' data.
% 
% NOTES: When using 'DMT' the process is a bit slow because each miniseed
% file has to be loaded and the start time read from the header. This is
% because the file name does not contain information about the timing of
% the data.
%
% A database file is written in the '.mat' format and used later for
% parellel processing of raw seismic data followed by correlation.
%
% INPUT:
%   project_directory = path to folder where output files will be written
%   data_directory    = path to data (not down into the NET/STA/ folders)
%   database_name     = name of the datebase file we create
%   file_type         = 'sac', 'seed', 'miniseed' (only 'sac' implemented)
%   data_structure    = 'SDS', 'BUD', 'IDDS', 'PDF', 'DMT' (only 'BUD' and 'DMT' currently implemented; 'DMT' only looks for data in the 'processed' folder!)
%   start_date        = first day of database ['YYYY-MM-DD HH:MM:SS.FFF']
%   end_date          = last day of database ['YYYY-MM-DD HH:MM:SS.FFF']
%   channel_list      = channel_list (cell list of channels to use; e.g. {'BHZ','BHE'})
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
% data_structure['DMT']  = "continuous???/NET.STA.LOC.CHAN"
% data_structure['ANT']  = "SAC/NET.STA.CHAN.YEARMONTHDAY"
%
% Written by: Dylan Mikesell (dylanmikesell@boisetate.edu)
% Last modified: 4 September 2017
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
% Check MATLAB version (currently working in MATLAB R2017a)
%--------------------------------------------------------------------------
disp( ['Using MATLAB release R' version('-release')] );
% assert( strcmp('2017a',version('-release'))==1,...
%     'Not using MATLAB 2017a. Current version of code works with 2017a.' );

assert( isdir( project_directory ) == 1,...
    'Project directory does not exist. Check directory PATH!')

%--------------------------------------------------------------------------
% Check data format
%--------------------------------------------------------------------------
assert(...
    strcmp(file_type,'sac')==1 ||...
    strcmp(file_type,'miniseed')==1 ||...
    strcmp(file_type,'seed')==1,...
    'Data type not recognized. Must be ''sac'', ''seed'' or ''miniseed''.');

if strcmp(file_type,'miniseed')
    file_type = 'seed'; % get correct name for waveform()
end

%--------------------------------------------------------------------------
% Check data_structure
%--------------------------------------------------------------------------
assert(...
    strcmp(data_structure,'BUD')==1 ||...
    strcmp(data_structure,'DMT')==1 ||...
    strcmp(data_structure,'ANT')==1,...
    'Data structure not recognized. Must be ''BUD'' or ''DMT'' or ''ANT''.');

%--------------------------------------------------------------------------
% scan through the data files in all subdirectories
%--------------------------------------------------------------------------

switch data_structure % account for different naming conventions
    
    case 'BUD' % data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
        % get list of files
        fileList = dir( fullfile( data_directory, '**/*.*') );
    case 'DMT' % data_structure['DMT'] = "continuous???/NET.STA.LOC.CHAN"
        % get list of files
        fileList = dir( fullfile( data_directory, '**/*.*') );
        % use data in 'processed' folder
        processed_idx = strfind( {fileList.folder}, 'processed' );
        kill_idx = cellfun('isempty', processed_idx);
        fileList(kill_idx) = [];
    case 'ANT' % data_structure['ANT'] = "SAC/NET.STA.CHAN.YEARMONTHDAY"
        % get list of files
        a = dir( fullfile( data_directory, '20*') );
        for this_dir = {a(1:numel(a)).name}
            fprintf('Scanning %s\n',fullfile( data_directory, this_dir{1}, 'SAC') );
            tmp_sac_files = dir( fullfile( data_directory, this_dir{1}, 'SAC', '*.sac') );
            fprintf('Found %d sac files.\n',numel(tmp_sac_files));
            if exist('fileList','var') ~= 1
                fileList = tmp_sac_files;
            else
                fileList = [fileList; tmp_sac_files];
            end
        end

end

% remove directories from the list
dirIdx = [fileList.isdir];
fileList(dirIdx) = [];
clear dirIdx;

% remove any hidden files that start with '.'
fileList = fileList( arrayfun( @(x) ~strcmp(x.name(1),'.'), fileList ) );

%--------------------------------------------------------------------------
% Determine unique networks and stations
%--------------------------------------------------------------------------
switch data_structure % account for different naming conventions
    
    case 'BUD' % data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
        
        % get information about the first file
        file_info      = split( fileList(1).name, '.' );
        file_station  = file_info{1};
        file_network  = file_info{2};
        file_location = file_info{3};
        file_channel  = file_info{4};

        % Build the unique station identifier
        instrument_list = {strcat(...
            char(file_network), '.', char(file_station), '.',...
            char(file_location), '.', char(file_channel) ) };
        num_instruments = 1;
        
        for iFile = 2 : numel( fileList )
            
            file_info     = split( fileList(iFile).name, '.' );
            file_station  = file_info{1};
            file_network  = file_info{2};
            file_location = file_info{3};
            file_channel  = file_info{4};

            % Build the unique station identifier
            instrument = strcat(...
                char(file_network), '.', char(file_station), '.',...
                char(file_location), '.', char(file_channel) );
            
            % Check is instrument is already in list
            if ~logical( sum( strcmp( instrument_list, instrument ) ) )
                fprintf('Adding %s\n', instrument);
                num_instruments = num_instruments + 1;
                instrument_list{ num_instruments } = instrument;
            end
            
        end % end loop over all file found
              
    case 'DMT' % data_structure['DMT'] = "continuous???/NET.STA.LOC.CHAN"
        
        % get information about the first file
        file_info     = split( fileList(1).name, '.' );
        file_network  = file_info{1};
        file_station  = file_info{2};
        file_location = file_info{3};
        file_channel  = file_info{4};
        
        % Build the unique station identifier
        instrument_list = {strcat(...
            char(file_network), '.', char(file_station), '.',...
            char(file_location), '.', char(file_channel) ) };
        num_instruments = 1;
        
        for iFile = 2 : numel( fileList )
            
            file_info     = split( fileList(iFile).name, '.' );
            file_network  = file_info{1};
            file_station  = file_info{2};
            file_location = file_info{3};
            file_channel  = file_info{4};
            
            % Build the unique station identifier
            instrument = strcat(...
                char(file_network), '.', char(file_station), '.',...
                char(file_location), '.', char(file_channel) );
            
            % Check is instrument is already in list
            if ~logical( sum( strcmp( instrument_list, instrument ) ) )
                fprintf('Adding %s\n', instrument);
                num_instruments = num_instruments + 1;
                instrument_list{ num_instruments } = instrument;
            end
        end % end loop over all files found 
        
    case 'ANT' % data_structure['ANT'] = "SAC/NET.STA.CHAN.YEARMONTHDAY"
                
        % get information about the first file
        file_info     = split( fileList(1).name, '.' );
        file_network  = file_info{1};
        file_station  = file_info{2};
        file_channel  = file_info{3};
        file_location = '00';
        
        % Build the unique station identifier
        instrument_list = {strcat(...
            char(file_network), '.', char(file_station), '.',...
            char(file_location), '.', char(file_channel) ) };
        num_instruments = 1;
        
        for iFile = 2 : numel( fileList )
            
            file_info     = split( fileList(iFile).name, '.' );
            file_network  = file_info{1};
            file_station  = file_info{2};
            file_channel  = file_info{3};
            file_location = '00';

            
            % Build the unique station identifier
            instrument = strcat(...
                char(file_network), '.', char(file_station), '.',...
                char(file_location), '.', char(file_channel) );
            
            % Check is instrument is already in list
            if ~logical( sum( strcmp( instrument_list, instrument ) ) )
                fprintf('Adding %s\n', instrument);
                num_instruments = num_instruments + 1;
                instrument_list{ num_instruments } = instrument;
            end
        end % end loop over all files found       
end

%--------------------------------------------------------------------------
% Check that channel is in channel_list
%--------------------------------------------------------------------------
disp('Removing channels not in channel list...')
tmp_split = split( instrument_list, '.' );
networks  = tmp_split(:,:,1);  
stations  = tmp_split(:,:,2); 
locations = tmp_split(:,:,3);
channels  = tmp_split(:,:,4);

keep_idx = zeros( numel(instrument_list), 1 ); % set all instruments to false
% Set instrument to true if channel matches channel list
for ii = 1 : numel( channel_list )
    idx = strcmp( channel_list{ii}, channels );
    keep_idx(idx) = 1;
end
% Keep only instruments that pass the channel test
instrument_list = instrument_list( logical( keep_idx ) );

%--------------------------------------------------------------------------
% Print information to the user
%--------------------------------------------------------------------------
tmp_split = split( instrument_list, '.' );
networks  = unique( tmp_split(:,:,1) );  
stations  = unique( tmp_split(:,:,2) ); 

fprintf('Found %d unique networks:\n', numel( networks ) );
fprintf('%s\n', networks{:} );
fprintf('Found %d unique stations:\n', numel( stations ) );
fprintf('%s\n', stations{:} );

%--------------------------------------------------------------------------
% Build a structure to keep track of days with data
%--------------------------------------------------------------------------
days      = datenum( start_date ) : datenum( end_date ); % integer list of days
nDays     = numel( days ); % number of days in the table
nStations = numel( instrument_list );

stationData                  = struct;
stationData.Date             = days;
stationData.DataTable        = cell( nStations, nDays );
stationData.DataTable(:)     = {'N'}; % 'N' stands for no data
stationData.files            = fileList;
stationData.dataDirectory    = data_directory;
stationData.projectDirectory = project_directory;
stationData.fileType         = file_type;

%--------------------------------------------------------------------------
% Populate the structure database with files
%--------------------------------------------------------------------------

data_count = 0; % keep track of how many data you actually add to table

for iFile = 1 : numel( fileList )
    
    fprintf( 'Processing: %s\n', fileList(iFile).name );
    
    switch data_structure % account for different naming conventions
        
        case 'BUD' % data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
            file_info       = split( fileList(iFile).name, '.' );
            file_station    = file_info{1};
            file_network    = file_info{2};
            file_location   = file_info{3};
            file_channel    = file_info{4};
            file_year       = file_info{5};
            file_julian_day = file_info{6};
        case 'DMT' % data_structure['DMT'] = "continuous???/NET.STA.LOC.CHAN"
            file_info     = split( fileList(iFile).name, '.' );
            file_network  = file_info{1};
            file_station  = file_info{2};
            file_location = file_info{3};
            file_channel  = file_info{4};    
        case 'ANT' % data_structure['DMT'] = "continuous???/NET.STA.LOC.CHAN"
            file_info     = split( fileList(iFile).name, '.' );
            file_network  = file_info{1};
            file_station  = file_info{2};
            file_channel  = file_info{3};
            file_location = '00';
            file_date     = file_info{4};
            file_year     = file_date(1:4);
            file_month    = file_date(5:6);
            file_day      = file_date(7:8);
            
    end

    % Build the unique station identifier
    instrument = strcat(...
        char(file_network), '.', char(file_station), '.',...
        char(file_location), '.', char(file_channel) );
    
    % check that instrument is in list
    stationIdx = strcmp( instrument, instrument_list );
    if sum( stationIdx ) == 1

        file_path = { fullfile( fileList(iFile).folder, fileList(iFile).name ) }; % get file with path
        
        % Have to actually read miniseed and seed data to get year and
        % julian day
        switch data_structure % account for different naming conventions
            case 'DMT'
                w = waveform( file_path{1}, 'seed' );
                % Keep in mind that when we download with obspyDMT, we add
                % 1 second to data file start
                dataDate = datenum( datestr(get(w,'start')+1/3600/24), 'dd-mmm-yyyy' );
            case 'BUD'
                % compute the data date
                dataDate = datenum( [file_year '.' file_julian_day], 'yyyy.dd' );
            case 'ANT'
                % compute the data date
                dataDate = datenum( [file_year '.' file_month '.' file_day], 'yyyy.mm.dd' );
                
        end
        
        % check that data is in time window of interest
        if dataDate <= datenum( end_date ) && dataDate >= datenum( start_date )

            dayIdx = dataDate - datenum( start_date ) + 1; % ordinal number
            
            stationIdx = find( stationIdx == 1); % determine instrument row in data table

            fprintf( 'Adding waveform %s on %s\n', ...
                instrument_list{stationIdx}, ...
                datestr( dayIdx-1 + datenum(start_date) ) );

            stationData.DataTable( stationIdx, dayIdx ) = file_path; % save file with path
            data_count = data_count + 1; % update data counter
        end % end if over date range
    end % end if over instrument list
end % end loop through files

%--------------------------------------------------------------------------
% Finish and save the database
%--------------------------------------------------------------------------
save( fullfile( project_directory, database_name ), 'stationData');
fprintf( 'Found %d data files.\n', numel( fileList ) );
fprintf( 'Added %d data files to table.\n', data_count );
fprintf( 'Finished writing database file: %s\n', fullfile( project_directory, database_name ) );

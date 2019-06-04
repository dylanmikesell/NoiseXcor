function stationData = initializeTable( project_directory, data_directory, ...
    database_name, file_type, data_structure, start_date, end_date, ...
    channel_list, coordinate_file )
%
% USAGE: stationData = initializeTable( project_directory, data_directory, ...
%   database_name, file_type, data_structure, start_date, end_date, ...
%   channel_list, coordinate_file )
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
%   coordinate_file   = txt or csv file containing station informaiton (1
%   row per station); each row should contain NAME, NET, LAT [deg], LON
%   [deg], ELE [m]
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

%--------------------------------------------------------------------------
% check inputs -- will raise errrors if bad inputs
%--------------------------------------------------------------------------
file_type = check_inputs( project_directory, data_directory, database_name,...
    file_type, data_structure, start_date, end_date, channel_list, ...
    coordinate_file);

%--------------------------------------------------------------------------
% scan through the data files in all subdirectories
%--------------------------------------------------------------------------
fileList = scan_files( data_structure, data_directory );

%--------------------------------------------------------------------------
% Determine unique networks and stations and get requested channels
%--------------------------------------------------------------------------
instrument_list = get_instruments( data_structure, fileList, channel_list );

% Window instrument list of testing/debugging
% instrument_list = instrument_list(1:6);
% Get just the ST* stations
% instrument_list = instrument_list(strncmp('YT.ST',instrument_list,5));

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
% Populate the structure database with files
%--------------------------------------------------------------------------
stationData = fill_database( fileList, data_structure, instrument_list, ...
    start_date, end_date, data_directory, project_directory, file_type, ...
    database_name, coordinate_file, channel_list );

end % initializeTable()
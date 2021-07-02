clear; close all; clc

% This example script scans the folder Example_SAC/DATA, which contains
% continuous seismic data from the TC and ZV networks surrounding Peteroa
% volcano in Chile/Argentina. The output of this script is a database
% called Peteroa_db.mat. This database is then load in the second script to
% compute correlations from all of the data in the found in DATA/.

% -------------------------------------------------------------------------
% The user needs to set the following variables.

% Must be complete path to where you want to save data
project_directory = '/Users/dmikesell/GIT/NoiseXcor/examples/Example_SAC';
% All files will be written to this directory.

% continuous data information (in the project folder for convenience)
data_directory = fullfile( project_directory, 'DATA'); % path to data
file_type       = 'sac'; % 'sac', 'seed', 'miniseed' (only 'sac' implemented)
data_structure  = 'BUD'; % 'SDS', 'BUD', 'IDDS', 'PDF', 'DMT' (only 'BUD', 'DMT', 'ANT' implemented)

% Station coodinates: this is a formatted text file
% Columns are:  NAME    NET     LAT     LON     ELE
coordinate_file = fullfile( data_directory, 'station_coordinates.csv'); 
% This is the coordinate file (CSV or TXT)

% enter time of earliest data -- MAT file names are relative to this.
start_date = '2012-01-01 00:00:00'; % ['YYYY-MM-DD HH:MM:SS.FFF']
end_date   = '2012-01-30 00:00:00'; % ['YYYY-MM-DD HH:MM:SS.FFF']

% The user can choose which channels are scanned for.
channel = {'BHZ','HHZ','EHZ'}; % a cell list of any channels to use

database_name  = 'Peteroa_db.mat'; % what do you want call your database

% -------------------------------------------------------------------------
% End of user input
% -------------------------------------------------------------------------

% create the data table
stationData = initializeTable( ...
    project_directory, ...
    data_directory, ...
    database_name, ...
    file_type, ...
    data_structure, ...
    start_date, ...
    end_date, ...
    channel, ...
    coordinate_file ...
    );
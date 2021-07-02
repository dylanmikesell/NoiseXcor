clear all
close all
clc

noisexPATH = '/Users/dmikesell/GIT/NoiseXcor/';
addpath( fullfile( noisexPATH, 'src' ) );

% This script will look into a single data directory and make a database of
% all possible data files that can be correlated.

% folder where all files will be written
project_directory = fullfile( noisexPATH, 'Example_SAC' );

data_directory = fullfile( project_directory, 'DATA'); % where are the data?

database_name  = 'Peteroa_db.mat'; % what do you want call your database

file_type       = 'sac'; % 'sac', 'seed', 'miniseed' (only 'sac' implemented)
data_structure  = 'BUD'; % 'SDS', 'BUD', 'IDDS', 'PDF', 'DMT' (only 'BUD', 'DMT', 'ANT' implemented)
channel         = {'BHZ','HHZ','EHZ'}; % a cell list of any channels to use

% Station coodinates: formatted text file
% Columns are:  NAME    NET     LAT     LON     ELE
coordinate_file = './station_coordinates.csv'; % This is the coordinate file (CSV or TXT)

% enter time of earliest data -- MAT file names are relative to this.
start_date = '2012-01-01 00:00:00'; % ['YYYY-MM-DD HH:MM:SS.FFF']
end_date   = '2012-01-30 00:00:00'; % ['YYYY-MM-DD HH:MM:SS.FFF']

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
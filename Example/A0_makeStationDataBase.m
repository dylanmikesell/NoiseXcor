clear all
close all
clc

% This script will look into a single data directory and make a database of
% all possible data files that can be correlated.

addpath('../src');

% folder where all files will be written
project_directory = '/hammer/SOFTWARE/NoiseXcor/Example';

data_directory = fullfile( project_directory, 'DATA'); % where are the data?

database_name  = 'Peteroa_db.mat'; % what do you want call your database

file_type       = 'sac'; % 'sac', 'seed', 'miniseed' (only 'sac' implemented)
data_structure = 'BUD'; % 'SDS', 'BUD', 'IDDS', 'PDF' (only 'BUD' implemented)

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
    end_date ...
    );
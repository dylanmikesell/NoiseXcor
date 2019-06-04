%--------------------------------------------------------------------------
% Populate the structure database with files
%--------------------------------------------------------------------------
function stationData = fill_database( fileList, data_structure, ...
    instrument_list, start_date, end_date, data_directory, ...
    project_directory, file_type, database_name, coordinate_file, ...
    channel_list)

% Build a structure to keep track of days with data
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
stationData.channel_list     = channel_list;

data_count = 0; % keep track of how many data you actually add to table

for iFile = 1 : numel( fileList )
    
%     fprintf( 'Processing: %s\n', fileList(iFile).name );
    
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
        case 'SDS' % data_structure['SDS']  = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
            file_info       = split( fileList(iFile).name, '.' );
            file_network    = file_info{1};
            file_station    = file_info{2};
            file_location   = '00';
            file_channel    = file_info{4};
            %file_type       = file_info{5};
            file_year       = file_info{6};
            file_julian_day = file_info{7};
            
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
            case 'SDS'
                % compute the data date
                dataDate = datenum( [file_year '.' file_julian_day], 'yyyy.dd' );
                
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
fprintf('Adding station location information.\n');
[latlon, stationName, network, elevation] = read_station_file( coordinate_file );

% Match the row of the data table with the station coordinate text file
tmp_split = split( instrument_list, '.' );
list_members = strcat( tmp_split(:,:,1), '.', tmp_split(:,:,2) );
coord_list = strcat( network, '.', stationName );
for i_inst = 1 : numel( instrument_list )
    geom_idx = strcmp( list_members(i_inst), coord_list);
    if ( sum(geom_idx)~= 0 )
        stationData.CoordTable( i_inst, 1 ) = latlon(geom_idx,1); % save lat
        stationData.CoordTable( i_inst, 2 ) = latlon(geom_idx,2); % save lon
        stationData.CoordTable( i_inst, 3 ) = elevation(geom_idx); % save elevation
        stationData.stationTag( i_inst ) =  list_members(i_inst); % save station name
    else % did not find the station coordinates
        stationData.CoordTable( i_inst, 1 ) = NaN; % save lat
        stationData.CoordTable( i_inst, 2 ) = NaN; % save lon
        stationData.CoordTable( i_inst, 3 ) = NaN; % save elevation
        stationData.stationTag( i_inst ) =  list_members(i_inst); % save station name
    end
end

% Finish and save the database
save( fullfile( project_directory, database_name ), 'stationData');
fprintf( 'Found %d data files.\n', numel( fileList ) );
fprintf( 'Added %d data files to table.\n', data_count );
fprintf( 'Finished writing database file: %s\n', fullfile( project_directory, database_name ) );

end % fill_database()
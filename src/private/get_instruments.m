%--------------------------------------------------------------------------
% Determine unique networks and stations
%--------------------------------------------------------------------------
function instrument_list = get_instruments( data_structure, fileList, ...
    channel_list)

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
             
    case 'SDS' % data_structure['SDS']  = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
        
        % get information about the first file
        file_info     = split( fileList(1).name, '.' );
        file_network  = file_info{1};
        file_station  = file_info{2};
        file_location = '00';
        file_channel  = file_info{4};
        file_type  = file_info{5};
        file_year  = file_info{6};
        file_jday  = file_info{7};
        
        % Build the unique station identifier
        instrument_list = {strcat(...
            char(file_network), '.', char(file_station), '.',...
            char(file_location), '.', char(file_channel) ) };
        num_instruments = 1;
        
        for iFile = 2 : numel( fileList )
            
            file_info     = split( fileList(iFile).name, '.' );
            file_network  = file_info{1};
            file_station  = file_info{2};
            file_location = '00';
            file_channel  = file_info{4};
            file_type  = file_info{5};
            file_year  = file_info{6};
            file_jday  = file_info{7};
            
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
        
end % switch

% Check that channel is in channel_list
disp('Removing channels not in channel list...');
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

end % get_instrument()
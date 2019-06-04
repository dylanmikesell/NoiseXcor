%--------------------------------------------------------------------------
% Load station coordinate information
%--------------------------------------------------------------------------
function [latlon, stationName, network, elevation] = readStationFile( filename )

[~,~,ext] = fileparts(filename);

switch lower(ext)
    case '.csv'
        
        fid = fopen( filename, 'r' );
        C   = textscan( fid, '%s%s%f%f%f','Delimiter',',');
        fclose(fid);
        
        stationName = C{1};
        network     = C{2};
        latlon      = [C{3}, C{4}];
        elevation   = C{5};
        
    case '.txt'
        fid = fopen( filename, 'r' );
        C   = textscan( fid, '%s%s%f%f%f');
        fclose(fid);
        
        stationName = C{1};
        network     = C{2};
        latlon      = [C{3}, C{4}];
        elevation   = C{5};       
end
function fileList = scan_files( data_structure, data_directory )

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
%         a = dir( fullfile( data_directory, '201002/*') );
        for this_dir = {a(1:numel(a)).name}
            fprintf('Scanning %s\n',fullfile( data_directory, this_dir{1}, 'SAC') );
            tmp_sac_files = dir( fullfile( data_directory, this_dir{1}, 'SAC', '*.sac') );
            fprintf('Found %d sac files.\n',numel(tmp_sac_files));
            if exist('fileList','var') ~= 1
                fileList = tmp_sac_files;
            else
                fileList = [fileList; tmp_sac_files];
            end % if
        end % for 
    case 'SDS' % "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
        % get list of files
        fileList = dir( fullfile( data_directory, '**/*.*') );

end % switch

% remove directories from the list
dirIdx = [fileList.isdir];
fileList(dirIdx) = [];
clear dirIdx;

% remove any hidden files that start with '.'
fileList = fileList( arrayfun( @(x) ~strcmp(x.name(1),'.'), fileList ) );

end % scan_files()

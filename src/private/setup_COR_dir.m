%--------------------------------------------------------------------------
% Setup the output directory
%--------------------------------------------------------------------------
function outputDir = setup_COR_dir( stationData, corrFilter )

cor_directory = fullfile( stationData.projectDirectory, 'COR' );

if exist(cor_directory,'dir') ~= 7 % directory does not exist
    [success,message,messageID] = mkdir(cor_directory); % make new directory
    if success
        fprintf('Created COR directory: %s\n', cor_directory);
    else
        disp(messageID);
        disp(message);
        error('Did not create COR directory. Check permissions');
    end
else % directory exists
    disp('Found COR directory in the project directory');
end

% Check whether or not we should overwrite the filter output directory
outputDir = fullfile( stationData.projectDirectory, ...
    'COR', num2str(corrFilter.filterNum, '%02d' ) );

if exist( outputDir, 'dir' ) ~= 7 % the directory does not exist
    [success,message,messageID] = mkdir(outputDir); % make new directory
    if success
        fprintf('Created filter directory: %s\n', outputDir);
    else
        disp(messageID);
        disp(message);
        error('Did not create filter directory. Check permissions');
    end    
else % the directory does exist
    
    % Determine if the user wants to overwrite the directory or not
    saveQuestion = questdlg(...
        ['Found an existing filter with number ', num2str(corrFilter.filterNum, '%02d' ),'. What do you want to do?'], ...
        'Filter folder already exists', ...
        'Overwrite','Append','Cancel','Cancel');
    
    % Handle response
    switch saveQuestion
        case 'Overwrite'
            disp(['Overwriting ' outputDir]);
            
            % Remove the directory and all its contents first
            [success, message] = rmdir(outputDir,'s'); % recursively remove old directory
            if success
                disp('Successfully removed old directory.');
            else
                disp(message);
                error('Could not remove old directory. Check permissions');
            end
            % Make the new directory
            [success,message,messageID] = mkdir(outputDir); 
            if success
                fprintf('Created filter directory: %s\n', outputDir);
            else
                disp(messageID);
                disp(message);
                error('Did not create filter directory. Check permissions');
            end
        case 'Cancel'
            disp('Canceling this job.');
            return % go back to MATLAB prompt after ending processing
        case 'Append'
            disp('Appending data. Caution -- existing data may be overwritten.')
    end
end

end % setup_COR_dir()
function [success,message,messageID] = checkOutputDir(outputDirectory)
%
% This script checks to see if the outputDirectory already exists. If
% it exists, the user is prompted to see whether or not they want to
% overwrite the directory.


% if success
%     fprintf('Created OUTPUT directory: %s\n', outputDir)
% else
%     
% end


directoryCreate = 0;

if exist( outputDirectory, 'dir' )
    
    % Determine if the user wants to overwrite the directory or not
    saveQuestion = questdlg(...
        'Do you want to overwrite the entire directory or just update the files?', ...
        'Directory already exists', ...
        'Overwrite','Update','Cancel','Update');
    
    % Handle response
    switch saveQuestion
        case 'Overwrite'
            disp(['Overwriting ' outputDirectory]);
            directoryCreate = 1;
        case 'Update'
            disp(['Will just update files in ' outputDirectory]);
        case 'Cancel'
            error('Canceling job');
    end
    
%     userInput = input(['Directory already exists - overwrite ' outputDirectory '? [Y/n]: '],'s');
%     if strcmp(userInput,'Y');
%         directoryCreate = 1;
%     elseif strcmp(userInput,'n');
%         disp('Not overwriting directory.');
%     else
%         userInput = input('Try again - enter [Y/n]');
%         if strcmp(userInput,'Y');
%             directoryCreate = 1;
%         elseif strcmp(userInput,'n');
%             disp('Not overwriting directory.');
%         else
%             disp('Learn how to enter yes or no!')
%             return
%         end
%     end
    
else % directory does not exist so just go ahead and create

    directoryCreate = 1;

end

if directoryCreate
    [success, message]          = rmdir(outputDirectory,'s'); % recursively remove old directory
    [success,message,messageID] = mkdir(outputDirectory); % make new directory
else
    success   = 0;
    message   = 0;
    messageID = 0;
end

return


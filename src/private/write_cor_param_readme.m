function write_cor_param_readme( outputDir, corrParam, corrFilter, ...
    dataBaseName )
% 
% USAGE: function writeCorParamReadme( outputDir, corrParam, corrFilter, ...
%     dataBaseName )
    
outputFile = fullfile( outputDir, 'README.txt' );

if ( exist( outputFile, 'file' ) == 2 )
    % file exists...so delete
    delete(outputFile);
end

fid = fopen( outputFile, 'w' );
fprintf( fid, 'dataBase Name: %s\n', dataBaseName );
fprintf( fid, 'filter No.:    %d\n', corrFilter.filterNum );
fprintf( fid, 'min frequency [Hz]:  %0.2f\n', corrFilter.fmin );
fprintf( fid, 'max frequency [Hz]:  %0.2f\n', corrFilter.fmax );
fprintf( fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n' );
fprintf( fid, '%% Sampling Params      %%\n' );
fprintf( fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n' );
fprintf( fid, 'frequency [Hz]:      %0.2f\n', corrParam.resampleFrequency );
fprintf( fid, 'corr length [min]:   %0.2f\n', corrParam.windowLengthMinutes );
fprintf( fid, 'overlap [%%]:         %0.2f\n', corrParam.overlapPercent );
fprintf( fid, 'output length [sec]: %0.2f\n', corrParam.tMaxOut );
fprintf( fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n' );
fprintf( fid, '%% Whitening Params     %%\n' );
fprintf( fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n' );
fprintf( fid, 'whiten flag:   %d\n', corrFilter.whiten );
fprintf( fid, 'whiten method: %s\n', corrFilter.wMethod  );
fprintf( fid, 'fmin [Hz]:     %0.2f\n', corrFilter.wfmin );
fprintf( fid, 'fmax [Hz]:     %0.2f\n', corrFilter.wfmax  );

% if not FTN then other amplitude normalization can occur
if ~strcmp(corrFilter.wMethod ,'ftn')
fprintf( fid, 'amp_norm flag: %d\n', corrFilter.ampNorm );
fprintf( fid, 'norm method]:  %s\n', corrFilter.timeNorm  );
end

fclose( fid );

end
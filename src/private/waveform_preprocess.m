function wCut = waveform_preprocess( windowStart, windowEnd, wCut, corrFilter, Fs, nSampWin, jdebug )
%----------------------------------------------------------------------
% Set default that whitening has NOT occurred
wCut(:) = addfield( wCut(:), 'isWhite', false );

% Apply spectral and amplitude normalization
if strcmp( corrFilter.wMethod, 'ftn' ) && corrFilter.whiten % apply frequency-time normalization in one step
    
    FB = [corrFilter.wfmin, corrFilter.wfmax]; % frequency band for whitening
    wCut  = waveform_whiten( wCut, FB, corrFilter.wMethod ); % whiten the waveforms
    
else % apply time normalization and frequency normalization
    
    % apply amplitude normalization
    if corrFilter.ampNorm
        wCut = waveform_normalization( wCut, corrFilter.timeNorm, corrFilter.wfmin );
    end
    
    if corrFilter.whiten
        FB = [corrFilter.wfmin, corrFilter.wfmax]; % frequency band for whitening
        wCut  = waveform_whiten( wCut, FB, corrFilter.wMethod ); % whiten the waveforms
    end
    
end
% pad traces to they start at correct time --> This must be done after
% pre-processing
wCut = pad( wCut, windowStart, windowEnd, 0 );

%----------------------------------------------------------------------
% pad waveforms to same length and align to start time
%     wCut = set( wCut, 'data_length', nSampWin );
%     wCut = align( wCut, datestr(windowStart,'dd/mm/yyyy HH:MM:SS.FFF'), Fs, 'pchip' );
%     wCut = set( wCut, 'data_length', nSampWin );
%     wCut = align( wCut, datestr(windowStart,'dd/mm/yyyy HH:MM:SS.FFF'), Fs );
%     plot(wCut); grid on;
%     xlim([10 20]);

%----------------------------------------------------------------------
% Interpolate all waveforms to align at start time of this window
timeVector = (0:nSampWin-1) / Fs; % the true time vector we want
traceStart = datestr( get(wCut,'start'), 'SS.FFF'); % get all start times

for iTrace = 1 : numel( wCut )
    
    if jdebug % print new start times afer alignment
        fprintf( 'start time: %s\n', traceStart(iTrace,:) );
    end
    
    % trace time vector depends on start time
    traceTimeVector = str2double(traceStart(iTrace,:)) + timeVector;
    data = double( wCut(iTrace) ); % get data
    data = interp1( traceTimeVector, data, timeVector, 'pchip', 0 ); % extrapolate with 0
    
    wCut(iTrace) = set( wCut(iTrace), 'data', data ); % update data
    wCut(iTrace) = set( wCut(iTrace), 'start', windowStart ); % update start time
    
    % and for good measure... update the waveform's history
    wCut(iTrace) = addhistory( wCut(iTrace), 'Aligned traces with 1D interpolation' );
end
%     plot(wCut); grid on;
%     xlim([10 20]);
%
end
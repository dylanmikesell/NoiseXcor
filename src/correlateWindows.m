function correlateWindows( W, corrParam, corrFilter, outputDirectory,...
    startDay)

% idx1 = strcmp(get(W,'station'),'CRIZ');  
% idx2 = strcmp(get(W,'station'),'PV01'); 
% idx = logical( idx1 + idx2 );
% W = W(idx);

jdebug = 0; % set to 1 for debug output

if isempty(W)
    disp('No data in this set of waveforms.');
    return
end

windowMin      = corrParam.windowLengthMinutes;
overlapPercent = corrParam.overlapPercent;
tMax           = corrParam.tMaxOut;
Fs             = corrParam.resampleFrequency;
% saveBeam       = corrParam.saveBeam;

%--------------------------------------------------------------------------
% setup window length
nptsDay  = 24 * 3600 * Fs + 1; % number of samples in 1 day
nSampWin = windowMin * 60 * Fs; % number of sample in the window
if overlapPercent == 0 % number of samples to move from one window to next
    nSlideWin = nSampWin;
else
    nSlideWin = floor( nSampWin * overlapPercent );
end
windowIdx  = 1 : nSlideWin : nptsDay; % starting index of windows
slideLimit = nptsDay - nSampWin; % check for window limit on back end
windowIdx( windowIdx > slideLimit ) = []; % remove windows that go over
nWindows   = numel( windowIdx ); % number of windows

fprintf('Number of windows %d.\n',nWindows);

%--------------------------------------------------------------------------
% start time information
% startTime = floor( min( get( W, 'start' ) ) );
startTime = startDay; % make sure you start on the day in the stationData table
startSTR = datestr(startTime,'YYYY-mm-DD HH:MM:SS.FFF');
fprintf('Correlations starting: %s\n', startSTR);

if jdebug % look at starting time of each waveform
    for thisWave = 1 : numel(W)
        startSTR2 = datestr(get( W(thisWave), 'start' ),'YYYY-mm-DD HH:MM:SS.FFF');
        fprintf('%s starts at %s\n', get( W(thisWave), 'station' ),startSTR2);
    end
end

%--------------------------------------------------------------------------
% setup output sampling
tMaxSamp = tMax * Fs; % number of samples to save
sampOut  = nSampWin-tMaxSamp+1 : nSampWin+tMaxSamp-1; % index for windowing correlations
nSampOut = numel(sampOut); % number of correlation samples to save

% loop over time windows in this day
for tt = 1 : nWindows
    
    windowStart = startTime + (windowIdx(tt)-1)/Fs/60/60/24;
    windowStartSTR = datestr(windowStart,'HH:MM:SS.FFF');
    
    windowEnd = startTime + (windowIdx(tt)-1+nSampWin)/Fs/60/60/24;
    windowEndSTR = datestr(windowEnd,'HH:MM:SS.FFF');
    
    fprintf('Correlating window %d of %d (%s: %s - %s)\n', tt, nWindows, datestr(startTime,'YYYY-mm-DD'), windowStartSTR, windowEndSTR );
    
    %----------------------------------------------------------------------
    % extract this time window
    wCut = extract( W, 'TIME', windowStart, windowEnd );
    % There will be zero samples if this window is outside the data range
%     wCut = align( wCut, datestr(windowStart,'dd/mm/yyyy HH:MM:SS.FFF'), get(wCut(1),'freq'), 'pchip' );
%     wCut = align( wCut, datestr(min(get(wCut,'start')),'dd/mm/yyyy HH:MM:SS.FFF'), get(wCut(1),'freq'), 'pchip' );
%     wCut = pad(wCut,windowStart,windowEnd,0);
    %----------------------------------------------------------------------
    % remove empty waveforms in this window and any that are all zeros
    % killIdx = ( get(wCut,'Data_Length') == 0 );
    killIdx = ( sum( double(wCut) ) == 0 ); % returns True is data_legnth=0 as well.
    wCut(killIdx) = []; % remove empty waveforms
    if isempty(wCut)
        disp('No data in this time window.');
        continue % go to next iteration of loop over windows
    end

    %----------------------------------------------------------------------
    % remove waveforms when duration is too short
    traceDuration = get( wCut, 'duration' ) .* (24*3600); % [s] duration
    freq = get(wCut,'freq');
    killIdx = ( ceil( freq .* traceDuration ) <= 2*floor( nSampWin * 0.05 ) ); 
    % returns True is duration is less than taper length
    wCut(killIdx) = []; % remove empty waveforms
    if isempty(wCut)
        disp('No data in this time window.');
        continue % go to next iteration of loop over windows
    end 

    %----------------------------------------------------------------------
    % process the remaining waveforms
    wCut = fillgaps( wCut, 0 ); % fill data gaps with zeros
    wCut = demean( wCut ); % remove mean
    wCut = detrend( wCut ); % remove linear trend

    %----------------------------------------------------------------------
    % Taper and resample waveforms before any other processing
    for iTrace = 1 : numel( wCut )
        % build 1 minute taper to make sure edges are zero
        nTap  = floor( nSampWin * 0.05 ); % taper 5% on each side of trace
        taper = tukeywin( 2*nTap-1, 1 ); % make taper that is twice as long
        
        data = double( wCut(iTrace) ); % get data
        data(1:nTap) = taper(1:nTap) .* data(1:nTap); % taper front
        nData = numel(data); % total number of samples
        data(nData-nTap+1:end) = taper(nTap:end) .* data(nData-nTap+1:end); % taper back
        
        resampleFactor = round( get( wCut(iTrace), 'freq' ) ) / Fs;
        remainder = mod( resampleFactor, floor(resampleFactor) );
        
        if remainder == 0
            % then resample as normal
            Q = resampleFactor;
            data = resample( data, 1, Q );  % see matlab's resample for specifics
        else % not an integer resample factor so first interpolate and then decimate
            multiplyFactor = round( 1 / remainder );
            data = resample( data, multiplyFactor, 1);
            Q = multiplyFactor * resampleFactor;
            data = resample( data, 1, Q );  % see matlab's resample for specifics
        end
        
        % put back into waveform, but don't forget to update the frequency too
        wCut(iTrace) = set( wCut(iTrace), 'data', data, 'freq', Fs );
        % and for good measure... update the waveform's history
        wCut(iTrace) = addhistory( wCut(iTrace), 'Resampled data outside of waveform object' );
        
    end

    %----------------------------------------------------------------------
    % Set default that whitening has NOT occurred
    wCut(:) = addfield( wCut(:), 'isWhite', false );
    
    % Apply spectral and amplitude normalization
    if strcmp( corrFilter.wMethod, 'ftn' ) && corrFilter.whiten % apply frequency-time normalization in one step
        
        FB = [corrFilter.wfmin, corrFilter.wfmax]; % frequency band for whitening
        wCut  = waveformWhiten( wCut, FB, corrFilter.wMethod ); % whiten the waveforms
        
    else % apply time normalization and frequency normalization
        
        % apply amplitude normalization
        if corrFilter.ampNorm
            wCut = waveformNormalization( wCut, corrFilter.timeNorm, corrFilter.wfmin );
        end
        
        if corrFilter.whiten
            FB = [corrFilter.wfmin, corrFilter.wfmax]; % frequency band for whitening
            wCut  = waveformWhiten( wCut, FB, corrFilter.wMethod ); % whiten the waveforms
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
%     if saveBeam
%         fprintf('Saving beamform data');
%     end
%     
    
    %----------------------------------------------------------------------
    % Double loop to cover all pairs of correlations
    nW = numel( wCut );
    for ii = 1 : nW
        
        WA = double( wCut(ii) ); % get data vector
        
        isWhitend = get( wCut(ii), 'isWhite' ); % check to see if data have been spectrally whitened already
        
        % loop over stations
        for jj = ii : nW
            
            if jdebug % print which stations being correlated
                fprintf('%s-%s\n', get(wCut(ii),'Station'), get(wCut(jj),'Station') );
            end
            
            WB = double( wCut(jj) ); % get data vector
            
            CC = normalizedCorrelation( WA, WB, Fs, corrParam.smoothMethod, corrParam.Wn, corrParam.K, isWhitend );
            % CC.c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
            % CC.c2: Simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
            % CC.c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
            
            % set the basic WAVEFORM properties
            statC = waveform(); % blank waveform object to store new correlations
            statC = set( statC, 'Station',  [get(W(ii),'Station') '-' get(W(jj),'Station')] );
            statC = set( statC, 'Channel',  [get(W(ii),'Channel') '-' get(W(jj),'Channel')] );
            statC = set( statC, 'Network',  [get(W(ii),'Network') '-' get(W(jj),'Network')] );
            statC = set( statC, 'Location', [get(W(ii),'Location') '-' get(W(jj),'Location')] );
            statC = set( statC, 'Start', windowStart );
            statC = set( statC, 'freq', Fs );
            statC = set( statC, 'Data_Length', nSampOut ) ;
            
            % add station location information
            % statC = addfield( statC, 'WALA', get( W(ii), 'STLA' ) );
            % statC = addfield( statC, 'WALO', get( W(ii), 'STLO' ) );
            % statC = addfield( statC, 'WAEL', get( W(ii), 'STEL' ) );
            % statC = addfield( statC, 'WBLA', get( W(jj), 'STLA' ) );
            % statC = addfield( statC, 'WBLO', get( W(jj), 'STLO' ) );
            % statC = addfield( statC, 'WBEL', get( W(jj), 'STEL' ) );
            
            % add windowed correlation functions and information
            if isfield( CC, 'c1' )
                statC = addfield( statC, 'c1', CC.c1(sampOut) );
            end
            if isfield( CC, 'c2' )
                statC = addfield( statC, 'c2', CC.c2(sampOut)  );
            end
            if isfield(CC,'c3')
                statC = addfield( statC, 'c3', CC.c3(sampOut)  );
            end
            statC = addfield( statC, 'smoothMethod', corrParam.smoothMethod );
            statC = addfield( statC, 'Wn', corrParam.Wn );
            statC = addfield( statC, 'K', corrParam.K );
            
            station   = get(statC,'station'); % get station pair
            startDate = [datestr(get(statC,'start'),'yyyy_mm_dd') '_window_' num2str(tt,'%03d')];
            
            if exist([outputDirectory '/' station],'dir') == 0
                mkdir([outputDirectory '/' station]);
            end
            
            fname = [outputDirectory '/' station '/' startDate '.mat'];
            save(fname,'statC','-v7.3'); % write NEW matrix file out
            
        end % jj loop over WB
    end % ii loop over WA
    
end % loop over time windows


end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
function Wout = waveformWhiten( W, F, wMethod )
%
% This is a wrapper function to whiten a single waveform trace based on
% Piero's BlanchMat.m routine that worked to whiten an entire matrix.
%
% USAGE: Wout = waveformWhiten( W, F )
%
% INPUT:
%   W = waveform object
%   F = frequency vector over which to whiten data [fmin,fmax] in Hz.
% Output:
%   Wout = whitened waveform object
%
% Written by Dylan Mikesell (dylanmikesell@boisestate.edu)
% Last modified 16 June 2016

nW   = numel( W ); % number of waveforms
Wout = W; % allocate waveform headers before updating the data field

% process each trace individually
for ii = 1 : nW
    
    %     fprintf( 'Whitening trace %d of %d\n', ii, nW );
    
    trace = double( W(ii) ); % get the data from waveform
    dt = 1 / get( W(ii), 'FREQ' ); % sample period of waveform
    
    % whiten data
    switch wMethod
        case 'ftn'
            trace = freqTimeNormalization( trace, F, dt );
        case 'poli'
            trace = whitenTrace( trace, F, dt );
        otherwise
            error('Only ftn and poli whitening approaches implement so far. Change wMethod')
    end
    
    %     [start_idx, end_idx] = findGaps( trace );
    %
    %     fprintf('Found %d data gap(s)\n', numel(start_idx) );
    %
    %     % loop through the non-zero sections of data
    %     for jj = 1 : numel( start_idx )
    %
    %         % get data slice
    %         tmp = trace( start_idx(jj) : end_idx(jj) );
    %
    %         % whiten data
    %         switch wMethod
    %             case 'ftn'
    %             trace( start_idx(jj) : end_idx(jj) ) = freqTimeNormalization( tmp, F, dt );
    %         case 'poli'
    %             trace( start_idx(jj) : end_idx(jj) ) = whitenTrace( tmp, F, dt );
    %         otherwise
    %             error('Only ftn and poli whitening approaches implement so far. Change wMethod')
    %         end
    %     end
    %
    % save the whitened waveform
    Wout(ii) = addhistory( Wout(ii), 'Whitened data over band [%0.2f, %0.2f] (Hz)', F(1), F(2) ); % add a history comment
    Wout(ii) = addhistory( Wout(ii), 'Whitened using %s method', wMethod ); % add a history comment
    Wout(ii) = set( Wout(ii), 'DATA', trace ); % replace old data with new
    Wout(ii) = addfield( Wout(ii), 'WhiteBand', F ); % add whitening band parameters
    Wout(ii) = addfield( Wout(ii), 'WhitenMethod', wMethod ); % add whitening band parameters
    Wout(ii) = set( Wout(ii), 'isWhite', true ); % indicate that whitening has occurred
    
end

end
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
function [start_idx, end_idx] = findGaps( trace )
%
% Usage: [start_idx, end_idx] = findGaps( trace )
%
% Description: Find data gaps that contains zeros
%
% Input:
%   trace = data vector
% Output:
%   start_idx = vector of start indices for data sections
%   end_idx = vector of end indices for data sections
%
% Example taken from here
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/13445
%
% Written by Dylan Mikesell (dylanmikesell@boisestate.edu)
% Last modified 7 July 2016

% data should be row vector
if size( trace, 2 ) == 1
    trace  = transpose( trace );
end

% identify gaps and process gaps separately
idx  = find( diff( [0 trace 0] ~= 0 ) );

% build 'good' data indices
start_idx = idx( 1 : 2 : end ); % take first index
end_idx   = idx( 2 : 2 : end ) - 1; % take second index

end
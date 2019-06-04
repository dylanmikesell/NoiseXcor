%--------------------------------------------------------------------------
% Cut, decimate/interpolate waveforms to time window
%--------------------------------------------------------------------------
function wCut = waveform_extract(windowStart,windowEnd, W, nSampWin, corrParam)

Fs = corrParam.resampleFrequency;

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
%     killIdx = ( sum( double(wCut) ) == 0 ); % returns True if data_legnth=0 as well.
%     wCut(killIdx) = []; % remove empty waveforms
%     if isempty(wCut)
%         disp('No data in this time window.');
%         continue % go to next iteration of loop over windows
%     end

%----------------------------------------------------------------------
% remove waveforms when duration is too short
traceDuration = get( wCut, 'duration' ) .* (24*3600); % [s] duration
freq = get(wCut,'freq');
killIdx = ( ceil( freq .* traceDuration ) <= 2*floor( nSampWin * 0.05 ) );
% returns True is duration is less than taper length
wCut(killIdx) = []; % remove empty waveforms
if isempty(wCut)
    disp('No data in this time window.');
    %         continue % go to next iteration of loop over windows
    return
end

%----------------------------------------------------------------------
% process the remaining waveforms in this time window
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
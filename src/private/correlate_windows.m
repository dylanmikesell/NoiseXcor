function correlate_windows( W, corrParam, corrFilter, outputDirectory,...
    startDay, CoordTable, stationTag, channel_list)

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

% Eliminate repeated stations and coordinates 
[unq_station,unq_idx] = unique( stationTag );
CoordTable2 = CoordTable(unq_idx,:);

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
        continue % go to next iteration of loop over windows
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
            wCut = waveformNormalization( wCut, corrFilter.timeNorm, corrFilter.wfmin );
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
    
    if corrParam.saveBeam
        fprintf('Saving beamform data');
        
        save_beam(CoordTable,stationTag,wCut,corrParam)
        
    end
    %
    %
    % clc
    %----------------------------------------------------------------------
    nW = numel( wCut ); % number of waveforms
    local_tags = get( wCut, 'ChannelTag' ); % station tags loaded today (e.g. 'YT.ST13..BHN')
    % Make complete list of stations --> Network.Station (e.g. 'YT.ST13')
    tmp_tag = cell(1,nW);
    for i_tag = 1 : nW
        tmp_tag{i_tag} = [local_tags(i_tag).network,'.',local_tags(i_tag).station];
    end
    % Get a unique list of stations without worrying about channels
    stations = unique(tmp_tag);
    nStats = numel( stations );

    % Do rotation for a station pair and then correlation
    for ii = 1 : nStats
        
        % isWhitend = get( wCut(ii), 'isWhite' ); % check to see if data have been spectrally whitened already
        isWhitend = 0;
        
        % Get waveforms
        idx = strcmp(stations{ii},tmp_tag);
        WA = wCut(idx);
        
        % Get coordinates
        idx = strcmp(stations{ii},unq_station);
        lat_wa = CoordTable2(idx,1);
        lon_wa = CoordTable2(idx,2);
        ele_wa = CoordTable2(idx,3);
        assert(numel(lat_wa)==1,'Did not find unique latitude for station 1');
        assert(numel(lon_wa)==1,'Did not find unique longitude for station 1');
        
        % Arrange the channels and extract waveforms
        if isa( get(WA,'channel'), 'char' ) % do not transpose
            wa_channels = get(WA,'channel');
        else % then transpose
            wa_channels = cell2mat( get(WA,'channel')' );
        end
        % wa_channels = cell2mat( get( WA, 'channel' )' );
        wa_components = wa_channels( :, 3 ); % get the last term
        
        e_idx = find( wa_components == 'E' );
        n_idx = find( wa_components == 'N' );
        z_idx = find( wa_components == 'Z' );
        
        % data_wa_enz = double(WA([e_idx,n_idx,z_idx]));
        data_wa_enz = zeros(nSampWin,3); % allocate and fill
        if ~isempty(e_idx)
            data_wa_enz(:,1) = double( WA(e_idx) );
        end
        if ~isempty(n_idx)
            data_wa_enz(:,2) = double( WA(n_idx) );
        end
        if ~isempty(z_idx)
            data_wa_enz(:,3) = double( WA(z_idx) );
        end
        
        % loop over stations
        for jj = ii : nStats
            
            if jj~= ii % cross correlation
                fprintf('Computing %s-%s correlation.\n',stations{ii},stations{jj});
                
                % Get channel tags for this station
                idx2 = strcmp(stations{jj},tmp_tag);
                WB = wCut(idx2); % extract waveforms for this station (e.g. E,N,Z)

                % Get coordinates
                idx2 = strcmp(stations{jj},unq_station);
                lat_wb = CoordTable2(idx2,1);
                lon_wb = CoordTable2(idx2,2);
                ele_wb = CoordTable2(idx2,3);
                assert(numel(lat_wb)==1,'Did not find unique latitude for station 2');
                assert(numel(lon_wb)==1,'Did not find unique longitude for station 2');

                % Arrange the channels and extract waveforms
                if isa( get(WB,'channel'), 'char' ) % do not transpose
                    wb_channels = get(WB,'channel');
                else % then transpose
                    wb_channels = cell2mat( get(WB,'channel')' );
                end
                wb_components = wb_channels(:,3); % get the last term
                % could do a check here to make sure that wb_components has
                % length = 3.
                
                e_idx = find(wb_components == 'E');
                n_idx = find(wb_components == 'N');
                z_idx = find(wb_components == 'Z');
                
                % data_wb_enz = double(WB([e_idx,n_idx,z_idx]));
                data_wb_enz = zeros(nSampWin,3); % allocate and fill
                if ~isempty(e_idx)
                    data_wb_enz(:,1) = double( WB(e_idx) );
                end
                if ~isempty(n_idx)
                    data_wb_enz(:,2) = double( WB(n_idx) );
                end
                if ~isempty(z_idx)
                    data_wb_enz(:,3) = double( WB(z_idx) );
                end

                
                % compute back-azimuth
                % BAZ = azimuth( 'rh', lat_wb, lon_wb, lat_wa, lon_wa, 'degrees' )
                % AZ = azimuth( 'rh', lat_wa, lon_wa, lat_wb, lon_wb, 'degrees' )
                BAZ = azimuth( 'gc', lat_wb, lon_wb, lat_wa, lon_wa, 'degrees' );
                % AZ = azimuth( 'gc', lat_wa, lon_wa, lat_wb, lon_wb, 'degrees' )
                
                data_wb_rtz = rotate_enz2rtz( data_wb_enz, BAZ, 'degrees' );
                data_wa_rtz = rotate_enz2rtz( data_wa_enz, BAZ, 'degrees' );
                
                % loop over correlation pairs requested by the user
                for i_comb = 1 : numel( corrParam.combinations )
                    fprintf('Computing %s correlation\n', corrParam.combinations{i_comb} );
                    
                    switch corrParam.combinations{i_comb}(1)
                        case 'R'
                            index1 = 1;
                        case 'T'
                            index1 = 2;
                        case 'Z'
                            index1 = 3;
                    end
                    
                    switch corrParam.combinations{i_comb}(2)
                        case 'R'
                            index2 = 1;
                        case 'T'
                            index2 = 2;
                        case 'Z'
                            index2 = 3;
                    end
                    
                    CC = normalizedCorrelation( data_wa_rtz(:,index1), data_wb_rtz(:,index2), Fs, corrParam.smoothMethod, corrParam.Wn, corrParam.K, isWhitend );
                    % CC.c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
                    % CC.c2: Simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
                    % CC.c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
                    
                    % Get receiver information
                    rec.name = get( WB(1),'station');
                    rec.chan = get( WB(1),'channel');
                    rec.chan(3) = corrParam.combinations{i_comb}(2);
                    rec.net = get( WB(1),'network');
                    rec.loc = get( WB(1),'location');
                    % WB = receiver = ST
                    
                    % Get source information
                    src.name = get( WA(1),'station');
                    src.chan = get( WA(1),'channel');
                    src.chan(3) = corrParam.combinations{i_comb}(1);
                    src.net = get( WA(1),'network');
                    src.loc = get( WA(1),'location');
                    % WA = source = EV
                    
                    
                    % set the basic WAVEFORM properties
                    statC = waveform(); % blank waveform object to store new correlations
                    
                    statC = set( statC, 'Start', windowStart );
                    statC = set( statC, 'freq', Fs );
                    statC = set( statC, 'Data_Length', nSampOut ) ;
                    statC = set( statC, 'Data', CC.c1(sampOut) );
                    
                    statC = set( statC, 'Station',  [src.name '-' rec.name] );
                    statC = set( statC, 'Channel',  [src.chan '-' rec.chan] );
                    statC = set( statC, 'Network',  [src.net '-' rec.net] );
                    statC = set( statC, 'Location', [src.loc '-' rec.loc] );
                    
                    
                    % add virtual source location information
                    statC = addfield( statC, 'EVLA', lat_wa );
                    statC = addfield( statC, 'EVLO', lon_wa );
                    statC = addfield( statC, 'EVEL', ele_wa );
                    % add receiver station location information
                    statC = addfield( statC, 'STLA', lat_wb );
                    statC = addfield( statC, 'STLO', lon_wb );
                    statC = addfield( statC, 'STEL', ele_wb);
                    % add inter-station information
                    statC = addfield( statC, 'BAZ', BAZ );
                    % DIST AZ (other SAC headers we could set)
                    
                    
                    %                 % add windowed correlation functions and information
                    %                 if isfield( CC, 'c1' )
                    %                     statC = addfield( statC, 'c1', CC.c1(sampOut) );
                    %                 end
                    %                 if isfield( CC, 'c2' )
                    %                     statC = addfield( statC, 'c2', CC.c2(sampOut)  );
                    %                 end
                    %                 if isfield(CC,'c3')
                    %                     statC = addfield( statC, 'c3', CC.c3(sampOut)  );
                    %                 end
                    %                 statC = addfield( statC, 'smoothMethod', corrParam.smoothMethod );
                    %                 statC = addfield( statC, 'Wn', corrParam.Wn );
                    %                 statC = addfield( statC, 'K', corrParam.K );
                    
                    station   = get(statC,'station'); % get station pair
                    startDate = [datestr(get(statC,'start'),'yyyy_mm_dd') '_window_' num2str(tt,'%03d')];
                    
                    basepath = fullfile(outputDirectory,corrParam.combinations{i_comb},station);
                    if ~isfolder(basepath)
                        mkdir(basepath);
                    end
                    
                    fname = fullfile(basepath,[startDate '.mat']);
                    save(fname,'statC','-v7.3'); % write NEW matrix file out
                    
                end % components ('ZZ','RR',etc.)
                
            else % jj=ii --> autocorrelation
                
                fprintf('Computing %s-%s autocorrelation.\n',stations{ii},stations{jj});
                
                % loop over autocorrelation pairs
                combinations = {'EE','NN','ZZ'};
                for i_comb = 1 : 3
                    
                    fprintf('Computing %s autocorrelation\n', combinations{i_comb} );
                    
                    CC = normalizedCorrelation( data_wa_enz(:,i_comb), data_wa_enz(:,i_comb), Fs, corrParam.smoothMethod, corrParam.Wn, corrParam.K, isWhitend );
                    % CC.c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
                    % CC.c2: Simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
                    % CC.c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
                    
%                     figure;
%                     plot(CC.c1)
                    
                    % Get receiver information
                    rec.name = get( WA(1),'station');
                    rec.chan = get( WA(1),'channel');
                    rec.chan(3) = combinations{i_comb}(2);
                    rec.net = get( WA(1),'network');
                    rec.loc = get( WA(1),'location');
                    % WB = receiver = ST
                    
                    % Get source information
                    src.name = get( WA(1),'station');
                    src.chan = get( WA(1),'channel');
                    src.chan(3) = combinations{i_comb}(1);
                    src.net = get( WA(1),'network');
                    src.loc = get( WA(1),'location');
                    % WA = source = EV
                    
                    % set the basic WAVEFORM properties
                    statC = waveform(); % blank waveform object to store new correlations
                    
                    statC = set( statC, 'Start', windowStart );
                    statC = set( statC, 'freq', Fs );
                    statC = set( statC, 'Data_Length', nSampOut ) ;
                    statC = set( statC, 'Data', CC.c1(sampOut) );
                    
                    statC = set( statC, 'Station',  [src.name '-' rec.name] );
                    statC = set( statC, 'Channel',  [src.chan '-' rec.chan] );
                    statC = set( statC, 'Network',  [src.net '-' rec.net] );
                    statC = set( statC, 'Location', [src.loc '-' rec.loc] );
                    
                    % add virtual source location information
                    statC = addfield( statC, 'EVLA', lat_wa );
                    statC = addfield( statC, 'EVLO', lon_wa );
                    statC = addfield( statC, 'EVEL', ele_wa );
                    % add receiver station location information
                    statC = addfield( statC, 'STLA', lat_wa );
                    statC = addfield( statC, 'STLO', lon_wa );
                    statC = addfield( statC, 'STEL', ele_wa);
                    % add inter-station information
                    statC = addfield( statC, 'BAZ', 0 );
                    % DIST AZ (other SAC headers we could set)
                    
                    station   = get(statC,'station'); % get station pair
                    startDate = [datestr(get(statC,'start'),'yyyy_mm_dd') '_window_' num2str(tt,'%03d')];
                    
                    basepath = fullfile(outputDirectory,combinations{i_comb},station);
                    if ~isfolder(basepath)
                        mkdir(basepath);
                    end
                    
                    fname = fullfile(basepath,[startDate '.mat']);
                    save(fname,'statC','-v7.3'); % write NEW matrix file out
                
                end % components ('ZZ','EE','NN')
                
            end % if autocorrelation or crosscorrelation
        end % jj loop over WB
    end %  ii loop over WA
end % loop over time windows
end % end correlateWindows() function

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------


% function save_beam( CoordTable, stationTag, wCut, corrParam )
%     %
%     % save_beam(CoordTable,stationTag,wCut)
%     %
%     % INPUT:
%     %   CoordTable = lat [deg], lon [deg], elevation [km]
%     %   stationTag = cell list of station names corresponds to row in
%     %   CoordTable --> YT.WISP
%     %   wCut = waveform object contain time series data
%
%     % Get data into frequency domain
%     dataFilt = double(wCut);
%     [npts, nstats] = size(dataFilt);
%
%     fs = get(wCut,'freq');
%     fs = unique(fs);
%     assert( numel(fs)==1, 'Waveforms do not have same dt! Cannot beamform.');
%     df   = fs/npts; % sample interval in fourier domain
%     fNyq = fs/2;    % Nyquist sampling frequency
%     DATA = fft(dataFilt, [], 1); % fourier transform each column of the matrix
%
%     % Get correct scaling for amplitudes
%     % Frequency vector creation and FFT scaling taken following:
%     % http://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/
%     DATA = DATA ./ npts; % scale the data by the number of points in FFT
%     % multiply all the non-unique amplitudes by 2
%     DATA( 2:ceil(npts/2), : ) = DATA( 2:ceil(npts/2), :) .* 2;
%
%     % take only the positive frequency part of the data
%     fArray = ( 0 : (npts-1) ) .* df; % frequency vector
%     fArray( fArray > fNyq ) = fArray( fArray > fNyq ) - ( fNyq * 2 );
%     fIdx = fArray >= 0; % get the indices of the positive frequencies
%     fPos = fArray( fIdx ); % frequency vector for positive frequencies
%     DATA = DATA(fIdx,:); % data for positive frequencies
%
%     % setup the slowness grid
%     cmin = corrParam.beam_cmin; % minimum phase
%     sxArray = linspace(-1/cmin,1/cmin,100);
%     syArray = linspace(-1/cmin,1/cmin,101);
%     [sx,sy] = meshgrid( sxArray, syArray ); % slowness grid
%
%     % get coordinates in local cartesian grid [km]
%     % lat0 = median(CoordTable(:,1))*(pi/180);
%     % lon0 = median(CoordTable(:,2))*(pi/180);
%     lat0 = mean(CoordTable(:,1));
%     lon0 = mean(CoordTable(:,2));
%     latlon = CoordTable(:,1:2);
%
%     % radius of the Earth
%     R = 6378.1; % [km]
%
%     % coordinates in the local cartesian system
%     for ii = 1 : size(latlon,1)
%         xyloc(ii,:) = (pi/180)*R*[ (latlon(ii,1)-lat0) (latlon(ii,2)-lon0)*cosd(lat0) ];
%         %     xyloc(ii,:) = R*[ ((latlon(ii,1)*(pi/180))-lat0) ((latlon(ii,2)*(pi/180))-lon0)*cos(lat0) ];
%     end
%
%     [~,freqIdx] = min(abs(fPos - 0.05));
%
%     w = 2 * pi * fPos(freqIdx); % [rad/s] omega at myFreq
%
%     ev = zeros( numel(syArray), numel(sxArray) ); % beam data
%     for ii = 1 : nstats % sum amplitudes of shifted stations
%         phi = ( sx .* xyloc(ii,1) ) + ( sy .* xyloc(ii,2) ); % phi = c_x*x + c_y*y
%         ev  = ev  + DATA(freqIdx,ii) * exp(1i * phi * w); % beamform
%     end
%     % compute beam power after normalizing by the number in the sum
%     ev = abs(ev ./ nstats).^2;
%     plot_beam(ev,cmin,sxArray,syArray,fPos(freqIdx))
%
%
%     % make example of array response at given slowness
%     sx0 = 0.0; % location in S_x
%     sy0 = -0.5; % location in S_y
%     [~,sxIdx] = min( abs( sxArray - sx0 ) ); % get the indices of synthetic source
%     [~,syIdx] = min( abs( syArray - sy0 ) );
%
%     ev = zeros( numel(syArray), numel(sxArray) ); % beam data
%     for ii = 1:nstats % sum amplitudes of shifted stations
%         phi = ( sx .* xyloc(ii,1) ) + ( sy .* xyloc(ii,2) ); % phi = c_x*x + c_y*y
%         phi0 = sxArray( sxIdx ) * xyloc(ii,1) + syArray( syIdx ) * xyloc(ii,2); % the initial phase from source at (sx0,sy0);
%         ev = ev + exp(1i * (phi+phi0) * w); % compute impulse response from signal with unit amplitude
%     end
%     % compute beam power after normalizing by the number in the sum
%     ev = abs(ev ./ nstats).^2;
%     plot_beam(ev,cmin,sxArray,syArray,fPos(freqIdx))
%
%
%     % compute the sensitivity matrix
%     sens_matrix = compute_sensitivity_matrix(sxArray,syArray,fPos(freqIdx),xyloc);
%
%
%
% end % save_beam

% function plot_beam(ev,cmin,sxArray,syArray,myFreq)
%
%     % Ray-parameter values
%     sMin = 0;
%     sMax = 1/cmin;
%     ds = (sMax-sMin) / 3; % we want three circles at different rayparameters
%     sArray = sMin : ds : sMax;
%
%     [sx,sy] = meshgrid( sxArray, syArray ); % slowness grid
%     s = sqrt( sx.^2 + sy.^2 ); % compute ray parameter at each node in grid
%     killIdx = (s >= sMax); % index of grid points with ray-parameter > rpMax
%     ev(killIdx) = NaN; % set these to zero and use yet_white colormap so they do not show up
%
%     h = figure('Color','w');
%     % plot the beamformed image
%     pcolor(syArray,sxArray,rot90(ev',2)); shading('interp'); hold on;
%     axis xy; axis('square'); axis(gca,'off'); % turn off x and y labels
%     axis(gca,sMax*[-1.15 1.15 -1.15 1.15]); % set axis limits
%     % colorbar and label
%     hh = colorbar;
%     set(get(hh,'YLabel'),'String','Beam Power'); % colorbar title
%
%     % plot radial spokes
%     th =  (0:5) * 2 * pi / 12; % we have 6 spokes
%     cst = cos(th);
%     snt = sin(th);
%     cs = [-cst; cst];
%     sn = [-snt; snt];
%     plot( sMax*cs, sMax*sn, 'k', 'LineWidth', 2 );
%
%     % annotate spokes in degrees
%     rt = 1.1 * sMax; % rayparameter max limit for annotations
%     for ii = 1 : length(th)
%         text( rt*cos(th(ii)), rt*sin(th(ii)), int2str(th(ii)*180/pi), 'HorizontalAlignment', 'Center' );
%         th2 = pi + th(ii); % shift by 180 to plot on other side as well
%         text( rt*cos(th2), rt*sin(th2), int2str(th2*180/pi), 'HorizontalAlignment', 'Center' );
%     end
%     text( 1.25*sMax, -0.25*sMax, 'Backazimuth [deg]' );
%     title(['Beamform output at ' num2str(myFreq) ' (Hz)'],'Position',[-1.2*sMax -1.2*sMax -1]);
%
%     % set view to clockwise from North
%     view([90,-90]);
%
%     % plot the circle for ray-parameter
%     dth = 1; % Azimuth values
%     thArray = ( 0 : dth : 360 - dth ) .* pi/180; % [rad]
%     for ii = 2 : numel(sArray)
%         plot( sArray(ii).*cos(thArray), sArray(ii).*sin(thArray), 'k', 'LineWidth', 2 );
%     end
%
%     %  annotate ray-parameters
%     cshift = cos( 45 * pi/180 ); % annotate at 45 degrees
%     sshift = sin( 45 * pi/180 );
%     for ii = 2 : numel(sArray)
%         text( (sArray(ii)+ds*0.05)*cshift, (sArray(ii)+ds*0.05)*sshift, [num2str(sArray(ii),'%2.2f')], 'VerticalAlignment', 'Bottom' );
%     end
%
% end % plot_beam
%
% function sens_matrix = compute_sensitivity_matrix(sxArray,syArray,myFreq,xyloc)
%
%     w = 2*pi*myFreq;
%     [sx,sy] = meshgrid( sxArray, syArray ); % slowness grid
%     nstats = size(xyloc,1);
%
%     nx = numel(sxArray);
%     ny = numel(syArray);
%     nxy = ny * ny; % total number of grid points
%     sens_matrix = zeros(nxy); % sensitivity matrix
%
%     for xx = 1 : nx % loop along X-direction
%         for yy = 1 : ny % loop along Y-direction
%
%             ev = zeros( numel(syArray), numel(sxArray) ); % make zero for this grid node
%             for ii = 1:nstats % sum amplitudes of shifted stations
%                 phi = ( sx .* xyloc(ii,1) ) + ( sy .* xyloc(ii,2) ); % phi = c_x*x + c_y*y
%                 phi0 = sxArray( xx ) * xyloc(ii,1) + syArray( yy ) * xyloc(ii,2); % the initial phase from source at (sx0,sy0);
%                 ev = ev + exp(1i * (phi+phi0) * w); % compute impulse response from signal with unit amplitude
%             end
%             % compute beam power after normalizing by the number in the sum
%             ev = abs(ev ./ nstats).^2;
%
%             % make into a vector and put in row of stacked power matrix
%             cnt = ( xx - 1 ) * ny + yy;
%             sens_matrix(cnt,:) = reshape( ev, 1, nxy ); % make a single row vector
%         end % yy
%     end % xx
%
% end % sens_matrix

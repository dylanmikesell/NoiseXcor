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
    windowEnd = startTime + (windowIdx(tt)-1+nSampWin)/Fs/60/60/24;
    fprintf('Correlating window %d of %d (%s: %s - %s)\n',...
        tt, nWindows, datestr(startTime,'YYYY-mm-DD'),...
        datestr(windowStart,'HH:MM:SS.FFF'),...
        datestr(windowEnd,'HH:MM:SS.FFF') );
    
    wCut = waveform_extract( windowStart, windowEnd, W, nSampWin, corrParam );
    
    wCut = waveform_preprocess( windowStart, windowEnd, wCut, corrFilter, Fs, nSampWin, jdebug );
    
    %     if corrParam.saveBeam
    %         fprintf('Saving beamform data');
    %
    %         save_beam(CoordTable,stationTag,wCut,corrParam)
    %
    %     end
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
                    
                    CC = normalized_correlation( data_wa_rtz(:,index1),...
                        data_wb_rtz(:,index2), Fs, corrParam.smoothMethod,...
                        corrParam.Wn, corrParam.K, isWhitend );
                    % CC.c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
                    % CC.c2: Simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
                    % CC.c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
                    
                    write_correlation( corrParam.combinations, i_comb, WB, WA, windowStart, Fs,...
                        nSampOut, CC, lat_wa, lon_wa, ele_wa, lat_wb, lon_wb, ele_wb, BAZ, ...
                        tt, outputDirectory, sampOut)
                    
                end % components ('ZZ','RR',etc.)
                
            else % jj=ii --> autocorrelation
                
                fprintf('Computing %s-%s autocorrelation.\n',stations{ii},stations{jj});
                
                % loop over autocorrelation pairs
                combinations = {'EE','NN','ZZ'};
                for i_comb = 1 : 3
                    
                    fprintf('Computing %s autocorrelation\n', combinations{i_comb} );
                    
                    CC = normalized_correlation( data_wa_enz(:,i_comb), data_wa_enz(:,i_comb), Fs, corrParam.smoothMethod, corrParam.Wn, corrParam.K, isWhitend );
                    % CC.c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
                    % CC.c2: Simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
                    % CC.c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
                    BAZ = 0;
                    write_correlation( combinations, i_comb, WA, WA, windowStart, Fs,...
                        nSampOut, CC, lat_wa, lon_wa, ele_wa, lat_wa, lon_wa, ele_wa, BAZ, ...
                        tt, outputDirectory, sampOut)
%                     
%                     %                     figure;
%                     %                     plot(CC.c1)
%                     
%                     % Get receiver information
%                     rec.name = get( WA(1),'station');
%                     rec.chan = get( WA(1),'channel');
%                     rec.chan(3) = combinations{i_comb}(2);
%                     rec.net = get( WA(1),'network');
%                     rec.loc = get( WA(1),'location');
%                     % WB = receiver = ST
%                     
%                     % Get source information
%                     src.name = get( WA(1),'station');
%                     src.chan = get( WA(1),'channel');
%                     src.chan(3) = combinations{i_comb}(1);
%                     src.net = get( WA(1),'network');
%                     src.loc = get( WA(1),'location');
%                     % WA = source = EV
%                     
%                     % set the basic WAVEFORM properties
%                     statC = waveform(); % blank waveform object to store new correlations
%                     
%                     statC = set( statC, 'Start', windowStart );
%                     statC = set( statC, 'freq', Fs );
%                     statC = set( statC, 'Data_Length', nSampOut ) ;
%                     statC = set( statC, 'Data', CC.c1(sampOut) );
%                     
%                     statC = set( statC, 'Station',  [src.name '-' rec.name] );
%                     statC = set( statC, 'Channel',  [src.chan '-' rec.chan] );
%                     statC = set( statC, 'Network',  [src.net '-' rec.net] );
%                     statC = set( statC, 'Location', [src.loc '-' rec.loc] );
%                     
%                     % add virtual source location information
%                     statC = addfield( statC, 'EVLA', lat_wa );
%                     statC = addfield( statC, 'EVLO', lon_wa );
%                     statC = addfield( statC, 'EVEL', ele_wa );
%                     % add receiver station location information
%                     statC = addfield( statC, 'STLA', lat_wa );
%                     statC = addfield( statC, 'STLO', lon_wa );
%                     statC = addfield( statC, 'STEL', ele_wa);
%                     % add inter-station information
%                     statC = addfield( statC, 'BAZ', 0 );
%                     % DIST AZ (other SAC headers we could set)
%                     
%                     station   = get(statC,'station'); % get station pair
%                     startDate = [datestr(get(statC,'start'),'yyyy_mm_dd') '_window_' num2str(tt,'%03d')];
%                     
%                     basepath = fullfile(outputDirectory,combinations{i_comb},station);
%                     if ~isfolder(basepath)
%                         mkdir(basepath);
%                     end
%                     
%                     fname = fullfile(basepath,[startDate '.mat']);
%                     save(fname,'statC','-v7.3'); % write NEW matrix file out
%                     
                end % components ('ZZ','EE','NN')
                
            end % if autocorrelation or crosscorrelation
        end % jj loop over WB
    end %  ii loop over WA
end % loop over time windows
end % end correlate_windows() function

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

function write_correlation( combinations, i_comb, WB, WA, windowStart, Fs,...
    nSampOut, CC, lat_wa, lon_wa, ele_wa, lat_wb, lon_wb, ele_wb, BAZ, ...
    tt, outputDirectory, sampOut)

% Get receiver information
rec.name    = get( WB(1),'station');
rec.chan    = get( WB(1),'channel');
rec.chan(3) = combinations{i_comb}(2);
rec.net     = get( WB(1),'network');
rec.loc     = get( WB(1),'location');
% WB = receiver = ST

% Get source information
src.name    = get( WA(1),'station');
src.chan    = get( WA(1),'channel');
src.chan(3) = combinations{i_comb}(1);
src.net     = get( WA(1),'network');
src.loc     = get( WA(1),'location');
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
% DIST, AZ (other SAC headers we could set)

% % add windowed correlation functions and information
% if isfield( CC, 'c1' )
%     statC = addfield( statC, 'c1', CC.c1(sampOut) );
% end
% if isfield( CC, 'c2' )
%     statC = addfield( statC, 'c2', CC.c2(sampOut)  );
% end
% if isfield(CC,'c3')
%     statC = addfield( statC, 'c3', CC.c3(sampOut)  );
% end
% statC = addfield( statC, 'smoothMethod', corrParam.smoothMethod );
% statC = addfield( statC, 'Wn', corrParam.Wn );
% statC = addfield( statC, 'K', corrParam.K );

station   = get(statC,'station'); % get station pair
startDate = [datestr(get(statC,'start'),'yyyy_mm_dd') '_window_' num2str(tt,'%03d')];

basepath = fullfile(outputDirectory,combinations{i_comb},station);
if ~isfolder(basepath)
    mkdir(basepath);
end

fname = fullfile(basepath,[startDate '.mat']);
save(fname,'statC','-v7.3'); % write NEW matrix file out

end % write_correlation() function
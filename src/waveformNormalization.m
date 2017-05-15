function Wout = waveformNormalization( W, normType, fmin )
%
% amplitude normalize waveforms
%
%   normType = 'abs', 'rms' or 'bit'

nW = numel(W); % number of waveforms
Wout = W; % allocate waveform headers before updating the data field

for ii = 1 : nW
    
    trace = double( W(ii) ); % get the data from waveform
    dt = 1 / get( W(ii), 'FREQ' ); % sample period of waveform

    % amplitude normalize data
    trace = amplitudeNormalization( trace, normType, dt, fmin ); % amplitude normalize
        
%     [start_idx, end_idx] = findGaps( trace );
%     
%     % loop through the non-zero sections of data
%     for jj = 1 : numel( start_idx )
%         
%         % get data slice
%         tmp = trace( start_idx(jj) : end_idx(jj) );
%         % amplitude normalize data
%         trace( start_idx(jj) : end_idx(jj) ) = amplitudeNormalization( tmp, normType, dt, fmin ); % amplitude normalize
%     end
    
    % save the normalized waveform
    Wout(ii) = addhistory( Wout(ii), 'Amplitude normalized using %s', normType ); % add a history comment
    Wout(ii) = set( Wout(ii), 'DATA', trace ); % replace old data with new
    
end

return
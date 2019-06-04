function Wout = waveform_whiten( W, F, wMethod )
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
            trace = frequency_time_normalization( trace, F, dt );
        case 'poli'
            trace = whiten_trace( trace, F, dt );
        otherwise
            error('Only ftn and poli whitening approaches implement so far. Change wMethod')
    end
    
    %     [start_idx, end_idx] = find_gaps( trace );
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


%--------------------------------------------------------------------------
% Find data gaps
%--------------------------------------------------------------------------
function [start_idx, end_idx] = find_gaps( trace )
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

end % find_gaps()
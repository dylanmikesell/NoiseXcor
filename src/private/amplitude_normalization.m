function d_out = amplitude_normalization( d_in, normType, dt, fmin )
%
% USAGE: d_out = amplitudeNormalization( d_in, normType, dt, fmin )
%
% Function normalizes the amplitudes over the length of a trace. The
% weight window length (for abs or rms) is chosen based on data frequency
% band --> fmin. 
%
% INPUT:
%   d_in = raw data vector amplitude normalization
%   normType = 'abs', 'rms' or 'bit'
%   dt   = [s] sample interval
%   fmin = [Hz] minimum frequency we are considering
%        
% OUTPUT:
%   d_out = amplitude normalized data
%

% Written by Dylan Mikesell (dylanmikesell@boisestate.edu)
% Last modified 15 June 2016
% Original example by Elmer Ruigrok, 2013

%--------------------------------------------------------------------------
% INPUT CHECK
%--------------------------------------------------------------------------

[ ns, ncol ] = size( d_in ); % ROW or COLUMN vector
isflip = 0; % a flag which says if vector is flipped
if (ns == 1) % then it is a column vector
    d_in   = transpose( d_in ); % make row vector
    isflip = 1; % set flag for flipping back at end
    ns     = ncol;
end

%--------------------------------------------------------------------------
% SETUP WINDOW LENGTH
%--------------------------------------------------------------------------
% The total fiter length is 2n+1 or less depending on where in the trace
% you are. See if statement below.

% In one direction the filter length in seconds equals 1/(4*fmin)
wl = 1 / 4 / fmin; % (s) window length
n  = ceil( wl / dt ); % (samp) window length

%--------------------------------------------------------------------------
% PROCESS TRACE
%--------------------------------------------------------------------------
switch normType
    case 'abs'

        dabs = abs( d_in ); % compute absolute value outside of loop for speed
        bw   = zeros( ns, 1 ); % allocate the weight vector for speed

        for k = 1 : ns % sliding window, one sample at a time
            
            % actual window length is 1/2/fmin when we go forward and backward
            kpn = k + n; % forward index
            kmn = k - n; % backward index
            
            % get the moving window indices
            if ( k > n && k <= ns - n ) % middle of trace
                ptIdx = kmn : kpn;
            elseif ( k <= n ) % front of trace
                ptIdx = 1 : kpn;
            elseif ( k > ns-n ) % end of trace
                ptIdx = kmn : ns;
            end
            
            % based on index, compute weight at k
            bw(k) = sum( dabs( ptIdx ) ) / numel( ptIdx );
        end
        d_out = d_in ./ bw; % normalize the trace by 'Bensen' factor
        
    case 'rms'
        
        dsqd = d_in.^ 2; % compute squared value outside of loop for speed
        bw   = zeros( ns, 1 ); % allocate the weight vector for speed

        for k = 1 : ns % sliding window, one sample at a time
            
            % actual window length is 1/2/fmin when we go forward and backward
            kpn = k + n; % forward index
            kmn = k - n; % backward index
            
            % get the moving window indices
            if ( k > n && k <= ns - n ) % middle of trace
                ptIdx = kmn : kpn;
            elseif (k <= n) % front of trace
                ptIdx = 1 : kpn;
            elseif (k > ns-n)
                ptIdx = kmn : ns; % end of trace
            end
            
            % based on index, compute weight at k
            bw(k) = sum( dsqd( ptIdx ) ) / numel( ptIdx );
        end
        d_out = d_in ./ sqrt( bw ); % normalize the trace by 'sqrt' factor
       
    case 'bit'
        
        d_out = sign( d_in ); % sign bit trace

    otherwise
        
        error('Unknown amplitude normalization type (choose: abs, rms or bit )');

end

d_out( isnan( d_out ) ) = 0; % check for any bw = 0 division causing NaN
d_out( isinf( d_out ) ) = 0; % check for any bw = 0 division causing inf

%--------------------------------------------------------------------------
% OUTPUT CHECK
%--------------------------------------------------------------------------

% flip to match original input trace
if ( isflip == 1 )
    d_out = transpose( d_out );
end

return
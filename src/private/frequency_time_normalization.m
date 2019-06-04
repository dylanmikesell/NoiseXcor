function ftnTrace = frequency_time_normalization( data, F, dt )

isflip = 0;
if size(data,1) == 1
    isflip = 1;
    data = transpose(data);
end

fl = F(1); % lowest frequency to whiten
fh = F(2); % highest frequency to whiten

df   = fl / 4; % The narrow frequency interval (df) is set to one-fourth of the lowest frequency

fkvec = fl : df : fh;   % frequency vector
nf    = numel( fkvec ); % nf is the number of narrow frequency bands
fnyq  = 1 / ( 2 * dt ); % [Hz] Nyquist frequency

% tmp = detrend( data, 'constant' ); % demean
% tmp = detrend( tmp, 'linear' ); % detrend
tmp = data;
npts = numel( tmp ); % data length

% taper based on the fh parameter
taperTime   = 1 / fh; % [s] length of the taper
taperSample = round( taperTime / dt ); % number of samples in the single sided taper
R           = taperSample / npts; % percentage of the taper relative to total time length
window      = tukeywin( npts, 2*R );
% tmp         = window .* tmp; % taper the data

btord = 2; % butterworth filter order (2 is recommended in the Shen et al. 2012 paper)

ftnTrace = zeros( npts, 1 ); % allocate

for ii = 1 : nf % loop through frequencies

    f1 = fl + df * ( ii - 1 );
    f2 = fl + df * ( ii );
    
    % build narrow band filter
    [ B, A ] = butter( btord, [f1/fnyq f2/fnyq], 'bandpass' );
    
    a = filtfilt( B, A, tmp ); % zero phase forward and backward filter
    
    env = sqrt( a.^2 + imag( hilbert(a) ).^2 ); % compute the envelope
    
    ftnTrace = ftnTrace + a ./ env; % sum each new filtered trace
    
end

ftnTrace = window .* ftnTrace; % make sure edges are zero

% make output the same shape as the input
if isflip
    ftnTrace = transpose(ftnTrace);
end

return
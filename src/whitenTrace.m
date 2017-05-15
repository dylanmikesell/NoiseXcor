function Awhite = whitenTrace(A, F, dt)
%
% This function spectrally whitens the trace A.
%
% USAGE: Awhite = whitenTrace(A, F, dt)
%
% INPUT: 
%   A  = the trace (npts,1)
%   F  = a vector of min and max frequencies [fmin,fmax]
%   dt = sample interval in seconds
% OUTPUT:
%   Awhite = the whitened trace (npts,1)
%
% Written by Piero Poli (Original was blanchmat.m)
% Modified by Dylan Mikesell (mikesell@mit.edu)
% Last modified 2 June 2014

npts = size(A,1);

isflip = 0;
if npts == 1
    isflip = 1;
    A = transpose(A);
end

%---------------------------------------------------------------------
% Original blanchmat.m function

df = 1/(npts*dt); % sample interval in frequency domain

% taper based on the F(1) parameter
taperTime   = 1 / F(1); % [s] length of the taper
taperSample = round( taperTime / dt ); % number of samples in the single sided taper
R           = taperSample / npts; % percentage of the taper relative to total time length
window      = tukeywin( npts, 2*R );
A           = window .* A; % taper the data

TF = fft(A); % FFT the trace

Nmin = floor( F(1)/df ) + 1; % sample nearest low end of spectrum
Nmax = floor( F(2)/df ) + 1; % sample nearest high end of spectrum

Napod = min( [20, floor( (Nmax-Nmin+1)/2 )] );

% get indices of positive frequencies for the whitening band
N1 = max( [1, Nmin-Napod] );
N2 = Nmin + Napod;
N3 = Nmax - Napod;
N4 = min( [floor((npts-1)/2), Nmax+Napod] );

% negative frequencies
N5 = max( [(npts - Nmax - Napod), floor((npts-1)/2)] );
N6 = npts - Nmax + Napod;
N7 = npts - Nmin - Napod;
N8 = min([npts - Nmin + Napod npts]);

Awhite = zeros( size(A) ); % allocate for speed

% whiten negative frequencies with cosine taper on edges
Awhite(N1:N2,:) = repmat(cos((pi/2:-pi/(2*(N2-N1)):0)).',1,size(A,2));
Awhite(N2:N3,:) = 1;
Awhite(N3:N4,:) = repmat(cos((0:pi/(2*(N4-N3)):pi/2)).',1,size(A,2));

% whiten positive frequencies with cosine taper on edges
Awhite(N5:N6,:) = repmat(cos((pi/2:-pi/(2*(N6-N5)):0)).',1,size(A,2));
Awhite(N6:N7,:) = 1;
Awhite(N7:N8,:) = repmat(cos((0:pi/(2*(N8-N7)):pi/2)).',1,size(A,2));

% Added by DM on Feb. 27, 2017 because sometimes abs(TF)=0 and blows up Awhite
zeroIdx = ( abs(TF) ~= 0 ); % find values that are not zero before division

Awhite(zeroIdx) = real( ifft( TF(zeroIdx) ./ abs(TF(zeroIdx)) .* Awhite(zeroIdx) ) ); % inverse FFT to time domain

Awhite = window .* Awhite; % make sure edges are zero

% make output the same shape as the input
if isflip
    Awhite = transpose(Awhite);
end

end

% From Matt Haney. We could implement this as well. Just need the whtnd.m
% fucntion. It is somewhere on Dylan's MATLAB GIT site.

%     % Applying whitening to the beam data for comparison with all-CC whitened (19/07/2014)
%     % SPECTRAL WHITENING
%     if whiten
%         for jj = 1 : nW
%             tmp = double( W(jj) );
%             tmp = whtnd( fmin/2, fmin, fmax, fmax*2, 1/Fs, tmp, fmin/4 );
%             W(jj)=set( W(jj), 'Data', tmp );
%         end
%     end

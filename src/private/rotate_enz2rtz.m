function rtz = rotate_enz2rtz( enz, AZ, units )
% 
% USAGE: rtz = rotate_enz2rtz( enz, AZ, units )
% 
% Rotates horizontal components of a seismogram.
% 
% Description: 
% The North- and East-Component of a seismogram will be rotated in Radial
% and Transversal Component. The angle is given as the back-azimuth, that is
% defined as the angle measured between the vector pointing from the station
% to the source and the vector pointing from the station to the North.
% 
% INPUT:
%   enz = data matrix (npts,3) ([e,n,z))
%   AZ = azimuth
%   units = units of AZ ('degrees' or 'radians')
% OUTPUT:
%   rtz = rotated data matrix (npts,3) ([r,t,z))
% 

% Example take from Obspy rotate_ne_rt()
% https://docs.obspy.org/_modules/obspy/signal/rotate.html#rotate_ne_rt

switch units
    case 'radians'
        AZ = AZ * 180/pi;
end

assert( AZ>=0 & AZ<360,...
    'Back Azimuth should be between 0 and 360 degrees.');
assert( numel( enz(:,1) ) == numel( enz(:,2) ),...
    'North and East components have different lengths.');

% copy data matrix
rtz = enz;
% compute new horizontal components
rtz(:,1) = - enz(:,1) * sind(AZ) - enz(:,2) * cosd(AZ);
rtz(:,2) = - enz(:,1) * cosd(AZ) + enz(:,2) * sind(AZ);

end
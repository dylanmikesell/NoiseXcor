clear; close all; clc;

% This script will plot a correlation gather and the stacked trace. The
% user needs to point to the station pair folder with the inputDir
% variable.

% inputDir = './COR/00/ZZ/TENZ-PV03';
inputDir = './COR/00/ZZ/CRIZ-TENZ';
% inputDir = './COR/00/ZZ/CRIZ-CRIZ';

time_limit = 30; % [s]

files = dir( inputDir );
idx = [files.isdir];
files(idx) = []; % remove any directories (i.e. '.' and '..')

% Load the correlation waveform files
for ii = 1 : numel( files )
    filename = fullfile( inputDir, files(ii).name );
    fprintf('Loading %s\n',filename);
    load( filename );
    Cmat(:,ii) = double( statC );
end

dt     = 1 / get( statC, 'freq' );
npts   = get( statC, 'data_length' );
time_array = ( (-npts+1) / 2 : (npts-1) / 2 ) .* dt;

% Plot the stacked trace and all of the individual correlations
figure;

subplot(1,7,1);
plot( mean(Cmat, 2), time_array );
ylim([-time_limit +time_limit]); grid on;
ylabel('Time (s)'); title('Mean');

subplot(1,7,2:7);
imagesc([],time_array,Cmat); axis xy;
ylim([-time_limit +time_limit]); colorbar;
xlabel('Correlation No.'); title('Correlation Gather');
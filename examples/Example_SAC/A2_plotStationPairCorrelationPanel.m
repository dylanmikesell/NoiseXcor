clear all
close all
clc

% plot a correlation panel
inputDir = './COR/00/ZZ/TENZ-PV03';
inputDir = './COR/00/ZZ/CRIZ-TENZ';
inputDir = './COR/00/ZZ/CRIZ-CRIZ';

time_limit = 30; % [s]

files = dir( inputDir );
idx = [files.isdir];
files(idx) = [];



for ii = 1 : numel( files )
    load( fullfile( inputDir, files(ii).name ) );
%     Cmat(:,ii) = get( statC, 'c1' );
    Cmat(:,ii) = double( statC );
end

dt     = 1 / get( statC, 'freq' );
npts   = get( statC, 'data_length' );
time_array = ( (-npts+1) / 2 : (npts-1) / 2 ) .* dt;

figure;
subplot(1,2,1);
imagesc([],time_array,Cmat);
ylim([-time_limit +time_limit]); colorbar;
subplot(1,2,2);
plot( time_array, mean(Cmat, 2) );
xlim([-time_limit +time_limit]); grid on;
% 
% %% compute some stats about the data
% 
% addpath('../src');
% 
% % fmin = 0.05; % [Hz]
% % fmax = 0.5; % [Hz]
% % btord = 2;
% 
% % corrMatrix = filter_dylan(Cmat, dt, fmin, fmax, btord );
% corrMatrix = Cmat;
% 
% rmsCorr = rms( corrMatrix, 1 );
% % rmsbar  = mean( rmsCorr );
% rmsbar  = median( rmsCorr ); % use median instead of mean
% % rmsstd  = std( rmsCorr );
% rmsstd  = mad( rmsCorr, 1); % use median absolute value instead of std. dev.
% 
% %% filter the data in a few different bands to look at convergence
% 
% 
% figure;
% subplot(1,2,1);
% imagesc(1:size(corrMatrix,2), tArray, corrMatrix);
% ylim([-timeLimit timeLimit]);
% subplot(1,2,2);
% plot(tArray, sum(corrMatrix,2));
% xlim([-timeLimit timeLimit]); grid on;
% 
% %% PWS approach
% 
% tic
% 
% corrMatrix = corrMatrix( 1:(npts+1)/2, : ); % take acausal side
% 
% % stack correlation matrix
% stack = sum(corrMatrix,2);
% 
% % S-transform of linear stack
% [ stranl, fvec ] = S_transform_FD_fullspec( stack, dt );
% 
% % make phase weight
% sumr = zeros(size(corrMatrix,1));
% 
% % number of correlations to use
% Nv = size(corrMatrix,2);
% % power of stacking as in Schimmel et al. (2010; GJI)
% pwr = 2;
% 
% f1 = 1i * 2 * pi * tArray(1:(npts+1)/2); % compute outside of loop for speed
% 
% for ii = 1 : Nv
%     
%     % S-transform of ii-th seismogram
%     [stran,fvec] = S_transform_FD_fullspec( corrMatrix( :, ii ), dt );
%     
%     % Equation 6 in Schimmel et al. (2010; GJI)
%     sumr = sumr + ( stran ./ abs(stran) ) .* exp(  fvec' * f1 );
%     ii
% end
% 
% % Equation 6 in Schimmel et al. (2010; GJI)
% sumr = abs( sumr / Nv ).^pwr;
% 
% % apply phase weight to linear stack to get phase-weighted-stack
% % Equation 7 in Schimmel et al. (2010; GJI)
% stranpws = sumr .* stranl;
% 
% % transform back to the time domain
% pws = S_transform_inverse_fullspec( stranpws, fvec );
% 
% toc
% 
% %%
% 
% figure;
% plot(tArray(1:(npts+1)/2), stack ,'r'); hold on;
% plot(tArray(1:(npts+1)/2), pws ,'k'); 
% xlim([-timeLimit timeLimit]); grid on;
% 
% % figure;
% % plot(tArray((npts+1)/2:end), stack ,'r'); hold on;
% % plot(tArray((npts+1)/2:end), pws ,'k'); 
% 
% %%
% 
% inputDir2 = '/hammer/DATA/Llaima/deconData/STACKS/01/REF/ZZ/';
% filename = 'LL_BVL_LL_TRL.SAC';
% 
% output = rsac( fullfile( inputDir2, filename ) );
% % msnoiseStack = output(:,2);
% 
% %%
% 
% fmin = 0.5; % [Hz]
% fmax = 1; % [Hz]
% btord = 2;
% 
% o1 = filter_dylan(output(:,2), output(2,1)-output(1,1), fmin, fmax, btord );
% 
% o2 = filter_dylan(pws, dt, fmin, fmax, btord );
% 
% o3 = filter_dylan(stack, dt, fmin, fmax, btord );
% 
% 
% h = figure;
% plot( output(:,1), o1./max(o1) ,'r'); hold on;
% plot(-tArray(1:(npts+1)/2), o2./max(o2) ,'k', 'Linewidth', 3); 
% plot(-tArray(1:(npts+1)/2), o3./max(o3) ,'b');
% legend('msnoise - linear','mikesell - pws', 'mikesell - linear'); legend boxoff;
% xlim([-10 timeLimit]); grid on;
% xlabel('Time [s]'); ylabel('Norm. Amp. [a.u.]');
% 
% set( h, 'Position', [100 100 2000 900] );
% set( h, 'PaperPositionMode', 'auto' );
% set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
% set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
% set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
% 
% print( h, '-dpng', 'MSNoiseComparison_0p5_1p0_Hz_fixedDecimate.png');
% 
% 
% %% 
% 
% clc
% 
% [latlon, stationName, component, elevation] = readStationFile( '/hammer/DATA/Llaima/stationData/LlaimaStations.csv' );
% 
% SourceStation = 'BVL';
% ReceiverStation = 'TRL';
% 
% % find source and receiver
% recIdx = strcmp( ReceiverStation, stationName );
% srcIdx = strcmp( SourceStation, stationName );
% % compute the interstation distance
% [arclen,az] = distance(...
%     latlon(srcIdx,1), latlon(srcIdx,2),...
%     latlon(recIdx,1), latlon(recIdx,2) );
% dist = deg2km(arclen,'earth'); % [km] interstation distance
% fprintf( 'Dist: %2.2f [km]\n', dist );
% 
% %%
% 
% vmin = 1500; % [m/s]
% vmax = 3500; % [m/s]
% dv = 10;
% alpha = 1;
% f1 = 0.05;
% f2 = 2;
% df = 0.01;
% 
% % [fvMatrix, fArray, vArray, traceSNR] = ftanTrace(flipud(stack),dt,f1,f2,df,vmin,vmax,dv,alpha,dist*1000,0);
% % [fvMatrix, fArray, vArray, traceSNR] = ftanTrace(flipud(pws),dt,f1,f2,df,vmin,vmax,dv,alpha,dist*1000,0);
% 
% [fvMatrix, fArray, vArray, traceSNR] = ftanTrace(output(12000:end,2), output(2,1)-output(1,1),f1,f2,df,vmin,vmax,dv,alpha,dist*1000,0);
% 
% 
% figure;
% imagesc( fArray, vArray, fvMatrix ); axis xy;
% 
% 
% 
% 
% 

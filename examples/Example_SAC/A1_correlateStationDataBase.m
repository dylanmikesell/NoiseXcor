clear; close all; clc

% -------------------------------------------------------------------------
% The user needs to set the following variables.

% Include full path if not in current directory.
dataBaseName  = './Peteroa_db.mat'; 

%--------------------------------------------------------------------------
corrFilter.filterNum = 0; % filter number is used to write correlations and keep track of completed jobs
corrFilter.fmin      = 0.01; % [Hz] low end of filter
corrFilter.fmax      = 0.99; % [Hz] high end of filter
corrFilter.prefilt   = 1; % 1 if prefilter, 0 if no prefilter
%--------------------------------------------------------------------------
% time domain normalization
corrFilter.ampNorm  = 1; % =1 if normalize, =0 if no amplitude normalization
corrFilter.timeNorm = 'abs'; % 'abs', 'bit', 'rms'
%--------------------------------------------------------------------------
% spectral whitening
corrFilter.whiten  = 1; % 0=off, 1=on
corrFilter.wfmin   = 0.05; % [Hz] low end of whitening
corrFilter.wfmax   = 1.0; % [Hz] high end of whitening
corrFilter.wMethod = 'ftn';
% wMethod = type of whitening to apply. Can be the following:
% 'ftn'   = frequency-time normalization (Shen et al. 2012).
% 'poli'  = moving window from Piero Poli
% 'sgn'   = unit amplitude for all frequencies as in msnoise (not yet implemented)
% 'haney' = smoothed moving window from Matt Haney (not yet implemented)
%--------------------------------------------------------------------------
% correlation parameters
corrParam.windowLengthMinutes = 60*4; % [min] correlation window length
corrParam.overlapPercent      = 0.0; % size of overlap [e.g. 1=100%, 0=0%, 0.5=50%]
corrParam.tMaxOut             = 120; % [sec] maximum time length to save of the correlation function (this is for one direction so actual correlation trace will be twice the this length)
corrParam.resampleFrequency   = 5; % [Hz] save computation time by resampling data prior to correlation...also saves disk space when writing output
corrParam.combinations        = {'ZZ','RR'}; % can be ZZ,RR,TT,ZR,ZT,RZ,TZ,RT,TR,EE,NN

% extra parameters for amplitude recovery of correlations (development
% code). Just leave as is for now.
corrParam.smoothMethod = 'taper'; % can be 'taper' or 'median'
corrParam.Wn           = 3;
corrParam.K            = 2*corrParam.Wn-1;

%--------------------------------------------------------------------------
% parallel computation parameters
np = 1; % number of processors to use

% -------------------------------------------------------------------------
% End of user input
% -------------------------------------------------------------------------

% Run the correlations
initializeCorrelation( ...
    dataBaseName, ...
    corrFilter, ...
    corrParam, ...
    np...
    );

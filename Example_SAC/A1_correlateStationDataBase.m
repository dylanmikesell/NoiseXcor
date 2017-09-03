clear all
close all
clc

addpath('../src'); % add source path with MATLAB functions

dataBaseName  = './Peteroa_db.mat'; % what did you want call your database. 
% Include full path if not in current directory.

%--------------------------------------------------------------------------
corrFilter.filterNum = 0; % filter number is used to write correlations and keep track of completed jobs
%--------------------------------------------------------------------------
% time domain normalization
corrFilter.ampNorm  = 1; % =1 if normalize, =0 if no amplitude normalization
corrFilter.timeNorm = 'bit'; % 'abs', 'bit', 'rms'
%--------------------------------------------------------------------------
% spectral whitening
corrFilter.whiten  = 1; % 0=off, 1=on
corrFilter.wfmin   = 0.5; % [Hz] low end of whitening
corrFilter.wfmax   = 5.0; % [Hz] high end of whitening
corrFilter.wMethod = 'ftn';
% wMethod = type of whitening to apply. Can be the following:
% 'ftn'   = frequency-time normalization (Shen et al. 2012).
% 'poli'  = moving window from Piero Poli
% 'sgn'   = unit amplitude for all frequencies as in msnoise (not yet implemented)
% 'haney' = smoothed moving window from Matt Haney (not yet implemented)
%--------------------------------------------------------------------------
% correlation parameters
corrParam.windowLengthMinutes = 60; % [min] correlation window length
corrParam.overlapPercent      = 0.0; % size of overlap [e.g. 1=100%, 0=0%, 0.5=50%]
corrParam.tMaxOut             = 120; % [sec] maximum time length to save of the correlation function (this is for one direction so actual correlation trace will be twice the this length)
corrParam.resampleFrequency   = 20; % [Hz] save computation time by resampling data prior to correlation...also saves disk space when writing output

% Beamforming data routines (not fully implemented yet)
corrParam.saveBeam            = 0; % 1=save, 0=do not save (writes beamform matrix)
corrParam.beam_cmin           = 1.0; % [km/s] minimum velocity for beam forming
corrParam.coordFile           = []; % Station Coordinate file

% extra parameters for amplitude recovery of correlations (development
% code). Just leave as is for now.
corrParam.smoothMethod = 'taper'; % can be 'taper' or 'median'
corrParam.Wn           = 3;
corrParam.K            = 2*corrParam.Wn-1;
%--------------------------------------------------------------------------
% parallel computation parameters
np = 1; % number of processors to use

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% END USER INPUT
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Run the correlations
initializeCorrelation( ...
    dataBaseName, ...
    corrFilter, ...
    corrParam, ...
    np...
    );

# How to run wavefield noise correlations with _NoiseXcor_

### The main codes:

The code works in two parts. 

* The first part (__A0_makeStationDataBase.m__) scans a folder for valid data and creates a database of files to be correlated.
* The second part (__A1_correlateStationDataBase.m__) reads the database, preprocesses the data and runs the correlations based on user specified parameters.

You can download test data [here](http://cgiss.boisestate.edu/~dmikesell/data.zip). Unzip the folder and put the data in the Example/ folder.

Each correlation is written to a new _waveform_ object with proper metadata. These correlation-derived waveforms can then be loaded into MATLAB for further processing (e.g. stacking, velocity analysis, etc.).

NOTE: There is a simple script (__A2_plotStationPairCorrelationPanel.m__) that can be used to plot a correlation panel.

---

### User parameters:

Within each of the main codes, the user needs to set a few parameters. Below we list each parameter and give some details. Note that a working example can be found in the Example/ folder.

1) __A0_makeStationDataBase.m__

	% choose a folder where all files will be written
	project_directory = '/home/NoiseXcor/Example';
	
	% indicate the data location
	data_directory = fullfile( project_directory, 'DATA');
	
	% give your database a name
	database_name  = 'Peteroa_db.mat'; 
	
	% Indicate the file format (e.g. 'sac', 'seed', 'miniseed')
	% This will eventually be any format that GISMOTOOLS can read
	file_type = 'sac'; % only 'sac' is implemented
	
	% we decided to follow the structure of the MSNoise data format
	% (e.g. ) 'SDS', 'BUD', 'IDDS', 'PDF'
	data_structure = 'BUD'; % only 'BUD' is currently implemented
	
	% enter time of earliest data to analyze
	start_date = '2012-01-01 00:00:00'; % ['YYYY-MM-DD HH:MM:SS.FFF']
	% enter time of latest data to analyze
	end_date   = '2012-01-30 00:00:00'; % ['YYYY-MM-DD HH:MM:SS.FFF']

2) __A1_correlateStationDataBase.m__

	%-----------------------------------------------------------------------
	
	% parallel computation parameters
	
	np = 1; % number of processors to use
	
	%-----------------------------------------------------------------------
	
	% filter number is used to keep track of completed jobs
	
	corrFilter.filterNum = 0; % can be any number
	
	%-----------------------------------------------------------------------
	
	% time domain normalization choice
	
	% =1 if normalize, =0 if no amplitude normalization
	
	corrFilter.ampNorm  = 1;
	 
	% 'abs', 'bit', 'rms'
	
	corrFilter.timeNorm = 'bit'; 
	
	%-----------------------------------------------------------------------
	
	% spectral whitening choice
	
	
	corrFilter.whiten  = 1; % 0=off, 1=on
	
	corrFilter.wfmin   = 0.5; % [Hz] low end of whitening
	
	corrFilter.wfmax   = 5.0; % [Hz] high end of whitening
	
	corrFilter.wMethod = 'ftn';
	
	% Notes on the type of whitening to apply (wMethod). This can be the followig:
	
	% 'ftn'   = frequency-time normalization (Shen et al. 2012).
	
	% 'poli'  = moving window from Piero Poli
	
	% 'sgn'   = unit amplitude for all frequencies as in msnoise (not yet implemented)
	
	% 'haney' = smoothed moving window from Matt Haney (not yet implemented)
	
	%-----------------------------------------------------------------------
	
	% correlation parameters choice
	
	% [min] correlation window length
	
	corrParam.windowLengthMinutes = 60;
	
	% size of overlap [e.g. 1=100%, 0=0%, 0.5=50%]
	
	corrParam.overlapPercent      = 0.0; 
	
	% [sec] maximum time length to save of the correlation function (this is for one direction so actual correlation trace will be twice the this length)
	
	corrParam.tMaxOut             = 120;
	
	% Save computation time by resampling data prior to correlation.
	
	% Also saves disk space when writing output waveforms.
	
	corrParam.resampleFrequency   = 20; % [Hz] 
	
	%-----------------------------------------------------------------------
	
	% beam forming parameters choice (not fully implemented yet)
	
	% 1=save, 0=do not save (writes beamform matrix)
	
	corrParam.saveBeam  = 0; 
	
	% [km/s] minimum velocity for beam forming
	
	corrParam.beam_cmin = 1.0;
	
	% [Hz] frequencies for beam forming, can be vector
	
	corrParam.beam_freqs = [];
	
	% Station Coordinate file 
	
	corrParam.coordFile = []; 
	
	% extra parameters for amplitude recovery of correlations (development code). Just leave as is for now.
	
	corrParam.smoothMethod = 'taper'; % can be 'taper' or 'median'
	
	corrParam.Wn           = 3;
	
	corrParam.K            = 2*corrParam.Wn-1;
	

3) __A2_plotStationPairCorrelationPanel.m__

	inputDir = './COR/00/CRIZ-TENZ';
	
You can give the input directory and the code will plot that station pair. Here we plot the correlations between station CRIZ and station TENZ.

### Extending to different formats

The data formats follow the MSNoise format.

	% data_structure['SDS']  = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
	% data_structure['BUD']  = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
	% data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
	% data_structure['PDF']  = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"

A user can alter src/initializeTable.m to incorporate new data formats. We have not yet had data in a format other than BUD, so we have not implemented anything else. If someone wants to share other formats, we would be happy to implement them. 

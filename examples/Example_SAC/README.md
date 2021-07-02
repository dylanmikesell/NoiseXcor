# NoiseXcor Example 1: continuous SAC data

This example demonstrates how to run wavefield noise correlations with SAC format data from two seismic networks in Chile.

### The main codes:

The code works in two parts. Each part has its own script.

1.  __A0_makeStationDataBase.m__ scans a folder for valid data and creates a database of files to be correlated.
2. __A1_correlateStationDataBase.m__ reads the database, preprocesses the data and runs the correlations based on user specified parameters.

You can download an example dataset [here](https://drive.google.com/file/d/1Wy4ycOEPMRJd8wzwEgbPuQ2T9rhkvu-J/view?usp=sharing). Unzip the file _sac_data.zip_. It should create a _DATA/_ folder. Put this _DATA/_ folder into the _Example\_SAC/_ directory.

### How does the code work?

Each correlation is stored as _waveform_ object (from GISMOTOOLS) with proper metadata (i.e. SAC header information). These correlation-derived waveforms can then be loaded into MATLAB for further processing (e.g. visualization, stacking, velocity analysis, etc.).

NOTE: There is a script (__A2_plotStationPairCorrelationPanel.m__) that can be used to plot a correlation panel (i.e. all of the correlations between two receivers).

---

### User parameters:

Within each of the main codes, the user needs to set a few parameters. Below we list each parameter and give some details.

1) __A0_makeStationDataBase.m__

This script sets up the inputs for the function initializeTable(). An example call is below.

```
initializeTable( ...
    project_directory, ...
    data_directory, ...
    database_name, ...
    file_type, ...
    data_structure, ...
    start_date, ...
    end_date, ...
    channel, ...
    coordinate_file ...
    );
```
Below is an explanation of the function inputs. The database can also be returned if a variable is assigned the output of initializeTable(). For example

```
stationData = initializeTable( ... );
```

### Function inputs

```
project_directory = full path to folder where output files will be written
data_directory    = full path to data
database_name     = name of the datebase file to be created
file_type         = 'sac', 'seed', 'miniseed' (only 'sac' implemented)
data_structure    = 'SDS', 'BUD', 'IDDS', 'PDF', 'DMT' (only 'BUD' and 'DMT' currently implemented; 'DMT' only looks for data in the 'processed' folder!)
start_date        = first day of database ['YYYY-MM-DD HH:MM:SS.FFF']
end_date          = last day of database ['YYYY-MM-DD HH:MM:SS.FFF']
channel_list      = channel_list (cell list of channels to use; e.g. {'BHZ','BHE'})
coordinate_file   = txt or csv file containing station informaiton (1 row per station); each row should contain NAME, NET, LAT [deg], LON [deg], ELE [m]
```

#### Coordinate file example (do NOT include a header line!)

```
PV03, ZV, -35.257, -70.502, 2448
AD2Z, TC, -35.148, -70.474, 2061
CRIZ, TC, -35.191, -70.524, 2909
TENZ, TC, -35.167, -70.504, 2264
```


### OUTPUT
```
%  stationData = a database structure containing file path information.
%   The database file is populated and also written to disk.
```


### Example

Here is the example from the _A0\_makeStationDataBase.m_ file in _Example\_SAC/_.

```
project_directory = '/Users/dmikesell/GIT/NoiseXcor/examples/Example_SAC';
data_directory = fullfile( project_directory, 'DATA');
file_type = 'sac';
data_structure  = 'BUD'; 
coordinate_file = fullfile( data_directory, 'station_coordinates.csv'); 
start_date = '2012-01-01 00:00:00';
end_date = '2012-01-30 00:00:00';
channel = {'BHZ','HHZ','EHZ'};
database_name  = 'Peteroa_db.mat';
```

2) __A1_correlateStationDataBase.m__

```
%-----------------------------------------------------------------------
	
% parallel computation
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
corrParam.resampleFrequency   = 20; % [Hz] 	
% Also saves disk space when writing output waveforms.
```	

#### A Note on the output folder structure

The individual waveforms are each written to .mat file that contains the _waveform_ object of the correlation. A waveform object has built in metadata that is useful later on for post processing. The files for each station pair are stored in folders for each station pair. The __COR__ directory is created in the _project_directory_ under that directory can be multiple "filter" directories.

3) __A2_plotStationPairCorrelationPanel.m__

```
inputDir = './COR/00/CRIZ-TENZ';
```
	
You can give the input directory and the code will plot that station pair. Here we plot the correlations between station CRIZ and station TENZ.



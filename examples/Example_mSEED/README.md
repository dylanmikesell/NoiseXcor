In this example, we use miniseed data downloaded from IRIS using the obspyDMT package. The data are stored in Example_mSEED/DATA/continuous1/processed. These data have been instrument deconvolve and are ready for correlation.

If you want to run the entire Polenet data set, use

	project_directory = '/Volumes/ESlab/Data_02/';
	data_directory = project_directory;
	
Otherwise you can change data_directory to a subset of the data with something like

	project_directory = '/Volumes/ESlab/Data_02/';
	data_directory = fullfile( project_directory, 'DATA2010/201002');


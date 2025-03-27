
// SERIES_SORTER: sorts the image and video captures from the experiment...
//only works if captures were saved with the naming scheme: OMTC_1, TFM_1, Cell_1, Null_1...


run("Bio-Formats Macro Extensions");

// USER INPUT 1/5: userDirectory requires the string of file directory name under which the experiemnt folder exists
userDirectory = "/D:/RAW DATA/"

// USER INPUT 2/5: experiment requires a string of the name of raw data folder from confocal
experiment = "ET_OMTCTFM_12kPa_SC_PC3CSS_PC31nM_48hour_mOrange_heat_13072022"

// USER INPUT 3/5: condition requires the name of folder for subexperiment 
condition = "PC3_1nMR1881_12kPa_SC"

directory_experiment = userDirectory + experiment

directory_condition = directory_experiment + "/" + condition

lif_file = directory_experiment +"/"+ experiment+ ".lif"

File.makeDirectory(directory_condition)

//USER INPUT 4/5: series_start requires the number corresponding to the FIRST series in condition
series_start = 11
//USER INPUT 5/5: series_start requires the number corresponding to the LAST series in condition
series_end = 20
//eg: in input examples above, positions 1-7 in the experiemnt were of LNCaP on 25 kPa
// series start to end corresponds to the position numbers and not the imagej series numbers 
//eg: for positions 1-5; even though imagej might display it as series 1-21, you would start to end as 1-5) 


for (series = series_start; series <= series_end; series++) {

	File.makeDirectory(directory_condition + "/series_" + series);
	save_tif (experiment,directory_condition, "OMTC", series);
	save_tif (experiment,directory_condition, "TFM", series);
	save_tif (experiment,directory_condition, "Cell", series);
	save_tif (experiment,directory_condition, "Null", series);
	save_tif (experiment,directory_condition, "AR", series);	
	

}

function save_tif (experiment, directory_condition, tif_type, series) {

	
	selectWindow(experiment + ".lif - " + tif_type + "_" + series);
		
	saveAs("Tiff", directory_condition + "/series_" + series + "/" + tif_type + "_" + series + ".tif");

	close(tif_type + "_" + series + ".tif" );
	
}	


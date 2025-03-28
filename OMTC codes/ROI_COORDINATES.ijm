//ROI_COORDINATES.ijm: finds coordinates for each ROI in a condition

//before running:
	//(1)Log window must be closed or empty
	//(2)any 512x512 image/video must be opened

//USER INPUT 1/4: experiment requires a string of the name of raw data folder from confocal
experiment = "EJ_OMTCTFM_LNCaP_12kPa_1nMR1881_CTGreen_48h_HEAT_06242022"

// USER INPUT 2/4: condition requires the name of folder for subexperiment 
condition = "LNCaP_1nMR1881_12kPa"

// USER INPUT 3/4: series_array requires an array of all the series to consider in the condition
series_array = newArray(1, 2 , 3, 4, 5);

// USER INPUT 4/4: userDirectory requires the string of file directory name under which the experiemnt folder exists
userDirectory = "/D:/RAW DATA/"

directory_experiment = userDirectory + experiment
directory_condition = directory_experiment + "/"+condition

for (h=1; h<=series_array.length; h++) {
	
	directory_series =directory_condition + "/series_" + series_array[h-1];
	directory_RoiSet = directory_series + "/RoiSet.zip";
	roiManager("Open", directory_RoiSet);
	
	roiCount = roiManager("count");
		for (i=1; i<=roiCount; i++) {
			roiManager("select", i-1);
			get_coordinates (i);
		
		}
	saveAs("Text", directory_series+"/RoiCoords.txt");
	roiManager("Deselect");
	roiManager("Delete");
	print("\\Clear");	

}


function get_coordinates (bead_num) {

	Roi.getCoordinates(x, y);
	
	  for (j=0; j<x.length; j++)
	    print("bead ",i, ",", x[j], ",", y[j]);
		
}



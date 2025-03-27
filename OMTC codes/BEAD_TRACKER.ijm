//BEAD_TRACKER: tracks lateral bead displacement 

//Instructions:
//1. USER INPUT (4): userDirectory, experiment, condition, series
//2. run script
//3. Trackmate 
//	a. input blob diameter (4.5 um) and threshold (5-12)manually 
// 	b. analyze spots, tracks, links
//	c. execute capture overlay
//	d. click OK
//	e. repeat for all ROIs

userDirectory ="/D:/RAW DATA/" //USER INPUT 1/4
experiment = "EJ_OMTCTFM_LNCaP_12kPa_1nMR1881_CTGreen_48h_HEAT_06242022" //USER INPUT 2/4
condition = "LNCaP_1nMR1881_12kPa" //USER INPUT 3/4
series = 6 //USER INPUT 4/4

directory_experiment = userDirectory + experiment
directory_condition = directory_experiment + "/"+condition
directory_series =directory_condition + "/series_" + series

roiCount = roiManager("count");
roi_array = newArray(0);

	for (i=1; i<=roiCount; i++) {

		bead_directory = directory_series+"/bead_"+ (i);
		File.makeDirectory(bead_directory);
		roiManager("select", i-1);  // select the ROI
        run("Duplicate...", "duplicate");
		run("Split Channels");
		
		//save C1 tiff
		selectWindow("C1-OMTC_"+ series + "-1.tif");
		C1_tiff = "C1_" + condition + "_series_" + series + "_bead_" + (i) + ".tif";
		print(C1_tiff);
		saveAs("Tiff", bead_directory + "/" + C1_tiff );
		close();
		
		run("Invert", "stack");
		run("Smooth", "stack");
		run("Enhance Contrast", "saturated=0.35");
		run("TrackMate");
		
		print(condition +"_series_" + series_num + "_bead_"+ bead_num+ ".csv");
		print(roicount);
		
		waitForUser("complete trackmate", "press OK once analysis is complete");
	

		selectWindow("TrackMate capture of C2-OMTC_" + series + "-1");
		saveAs("Tiff", bead_directory + "/C2_"+ condition +"_series_" + series + "_bead_" + (i) +".tif");
		close();
		close("C2_"+condition+"_OMTC_"+series+"_1.tif");
		close("C2-OMTC_"+series+"-1.tif");

		rename_roi (bead_directory, condition, series , i);

		selectWindow ("OMTC_" + series + ".tif");
}

Array.print(roi_array);
roiManager("Select", roi_array);
roiManager("Save", directory_series +"/RoiSet.zip");

close("*");
// saves ROI overlay 
function rename_roi (bead_directory, condition, series_num, bead_num) {

	bead_name = condition +"_series_" + series_num + "_bead_"+ bead_num;
	
	roiManager("Rename", bead_name );

	roi_array= Array.concat(roi_array, series_num-1);
	
	
	
}




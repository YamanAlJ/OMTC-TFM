function rename_roi (bead_directory, condition, series_num, bead_num) {

	bead_name = condition +"_series_" + series_num + "_bead_"+ bead_num;
	
	roiManager("Rename", bead_name );

	roi_array= Array.concat(roi_array, series_num-1);
userDirectory = "D:/RAW DATA/"

experiment = "EH_OMTCTFM_PC3AR_3kPa_10nMR1881_mOrange_48h_HEAT_06102022" 

condition = "PC3AR_10nMR1881_3kPa"

directory_experiment = userDirectory + experiment;
directory_condition = directory_experiment + "/"+condition;


series_array = newArray(2,3,4,5,6,7,8,9);
for (h=1; h<=series_array.length; h++) {
	
	series = series_array[h-1];
	directory_series =directory_condition + "/series_" + series;
	
	open(directory_series+"/TFM_"+ series+".tif");
	open(directory_series+"/Null_"+ series+".tif");
	open(directory_series+"/Cell_"+ series+".tif");
	
	
	position_dir = directory_condition + "/Series" + series; 
	File.makeDirectory(position_dir);
	
	selectWindow("Cell_"+ series+".tif"); 
	saveAs("Tiff", position_dir + "/Cell.tif"); 
	
	selectWindow("TFM_" + series+".tif");
	//run("Split Channels");
	//selectWindow("C1-" + "TFM_" + series+".tif");
	rename("Placeholder_Force.tif");								
	selectWindow("Null_" + series+".tif");
	//run("Split Channels");
	//selectWindow("C1-" + "Null_" + series+".tif");
	rename("Placeholder_Null.tif");

	
					
	// RUN BLEACH CORRECTOIN 
	run("Concatenate...", "  image1=Placeholder_Force.tif image2=Placeholder_Null.tif image3=[-- None --]");
	run("Properties...", "channels=1 slices=1 frames=2");
	rename("Bleach_Placeholder");
	 
	Stack.setChannel(1);
	rename("Bleach_placehold");
	
	run("Bleach Correction", "correction=[Histogram Matching]");
	rename("Merged");
	
	run("Make Substack...", "channels=1 frames=1"); //CHANGE with channel number
	rename("Placeholder_Force.tif");
	
	selectWindow("Merged");
	run("Make Substack...", "channels=1 frames=2"); // CHANGE
	rename("Null_Placeholder");
	
	run("Concatenate...", "  image1=Null_Placeholder image2=Placeholder_Force.tif image3=[-- None --]");
	close("\\Others");
	rename("MergedNull");
	
			
	// DEDRIFT 
	Stack.setChannel(1);
	Stack.setFrame(1)
	run("Properties...", "channels=1 slices=1 frames=2");
	run("Correct 3D drift", "channel=1 multi_time_scale sub_pixel edge_enhance only=0 lowest=1 highest=1 max_shift_x=20 max_shift_y=20 max_shift_z=0");
	saveAs("TIFF", position_dir + "/DeDrift_uncropped.tif");
	
	makeRectangle(0,0, 512, 512); // 1024
	run("Crop");
	
	makeRectangle(20,20, 472, 472); // 1024
	run("Crop");
	saveAs("TIFF", position_dir + "/NullAppended_DeDrift.tif");
	rename("NullAppended_DeDrift");
	
	// Step 5: Remove the null force image and reduce into a single image. 
	selectWindow("NullAppended_DeDrift"); 
	run("Make Substack...", "channels=1 frames=1");
	rename("Null Image");
	saveAs("Tiff", position_dir + "/Null.tif"); 
	close();
	
	File.makeDirectory(position_dir + "/Beads");
	selectWindow("NullAppended_DeDrift");
	run("Make Substack...", "channels=1 frames=2");
	rename("Beads_Forced");
	
	
	imagenum = 1;
	saveAs(position_dir + "/Beads/Beads_Forced_" + IJ.pad(imagenum, 4) +".tif");
						
	close("*");
	
}

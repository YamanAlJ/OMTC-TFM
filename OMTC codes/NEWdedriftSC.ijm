userDirectory = "D:/RAW DATA/"

experiment = "ET_OMTCTFM_12kPa_SC_PC3CSS_PC31nM_48hour_mOrange_heat_13072022" 

condition = "PC3_1nMR1881_12kPa_SC"

directory_experiment = userDirectory + experiment;
directory_condition = directory_experiment + "/"+condition;


series_array = newArray(11,12,13,14,15,16,17,18,19,20);
for (h=1; h<=series_array.length; h++) {
	
	series = series_array[h-1];
	directory_series =directory_condition + "/series_" + series;
	
	open(directory_series+"/TFM_"+ series+".tif");
	open(directory_series+"/Null_"+ series+".tif");
	open(directory_series+"/Cell_"+ series+".tif");
	open(directory_series+"/AR_"+ series+".tif");
	
	position_dir = directory_condition + "/Series" + series; 
	File.makeDirectory(position_dir);
	
	//selectWindow("Cell_"+ series+".tif"); 
	//saveAs("Tiff", position_dir + "/Cell.tif"); 
	
	selectWindow("TFM_" + series+".tif");
	//run("Split Channels");
	//selectWindow("C1-" + "TFM_" + series+".tif");
	rename("Placeholder_Force.tif");								
	selectWindow("Null_" + series+".tif");
	//run("Split Channels");
	//selectWindow("C1-" + "Null_" + series+".tif");
	rename("Placeholder_Null.tif");

	for (channum =1; channum <=3; channum++){ 
	selectWindow("Placeholder_Null.tif");
	run("Duplicate...", " title=" + channum);
	}
	
	//run("Merge Channels...", "c1=1 c2=2 c3=3 create ignore keep");		
	//run("Properties...", "channels=1 slices=3 frames=1");
	//rename("MergedNull");
	AR = "AR_" + series + ".tif";
	Cell = "Cell_"+series+".tif";
	run("Concatenate...", "image1=Placeholder_Force.tif image2=1 image3=[-- None --]");
	rename("Beads");
	run("Concatenate...", "image1=" + AR +  " image2=2 image3=[-- None --]");	
	rename("AR");
	run("Concatenate...", "image1=" + Cell +  " image2=3 image3=[-- None --]");	
	rename("Cell");	
	//run("Properties...", "channels=1 slices=3 frames=1");		
	//rename("Placeholder_Force");
	//run("Merge Channels...", "c1=Placeholder_Force c2=MergedNull create ignore keep");
	//rename("PlaceForceMergedNull");
	
					
	// RUN BLEACH CORRECTOIN 
	selectWindow("Beads");
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
	rename("MergedNullBEADS");
				//
	selectWindow("AR");
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
	rename("MergedNullAR");
				//
	selectWindow("Cell");
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
	rename("MergedNullCELL");
	
			
	// DEDRIFT 
	selectWindow("MergedNullBEADS");
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
	//
	selectWindow("MergedNullAR");
	Stack.setChannel(1);
	Stack.setFrame(1)
	run("Properties...", "channels=1 slices=1 frames=2");
	run("Correct 3D drift", "channel=1 multi_time_scale sub_pixel edge_enhance only=0 lowest=1 highest=1 max_shift_x=20 max_shift_y=20 max_shift_z=0");
	
	makeRectangle(0,0, 512, 512); // 1024
	run("Crop");
	
	makeRectangle(20,20, 472, 472); // 1024
	run("Crop");
	rename("NullAppended_DeDrift_AR");
	//
	selectWindow("MergedNullCELL");
	Stack.setChannel(1);
	Stack.setFrame(1)
	run("Properties...", "channels=1 slices=1 frames=2");
	run("Correct 3D drift", "channel=1 multi_time_scale sub_pixel edge_enhance only=0 lowest=1 highest=1 max_shift_x=20 max_shift_y=20 max_shift_z=0");
	
	makeRectangle(0,0, 512, 512); // 1024
	run("Crop");
	
	makeRectangle(20,20, 472, 472); // 1024
	run("Crop");
	rename("NullAppended_DeDrift_Cell");
	
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
	
	File.makeDirectory(position_dir + "/CellFluo");
	selectWindow("NullAppended_DeDrift_Cell");
	run("Make Substack...", "channels=1 frames=2");
	rename("CellFluo_Forced");
	
	imagenum = 1;
	saveAs(position_dir + "/CellFluo/CellFluo_Forced_" + IJ.pad(imagenum, 4) +".tif");
	
	File.makeDirectory(position_dir + "/AR");
	selectWindow("NullAppended_DeDrift_AR");
	run("Make Substack...", "channels=1 frames=2");
	rename("AR_Forced");
	
	imagenum = 1;
	saveAs(position_dir + "/AR/AR_Forced_" + IJ.pad(imagenum, 4) +".tif");
	
	
	
						
	close("*");
	
}

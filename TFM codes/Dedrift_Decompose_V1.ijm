
experiment = "EJ_OMTCTFM_LNCaP_12kPa_1nMR1881_CTGreen_48h_HEAT_06242022.lif" 
series = 1;

//Initialize 
selectWindow(experiment+" - TFM_" + series);
//run("Split Channels");
//selectWindow("C1-" + "TFM_" + series+".tif");
rename("Placeholder_Force.tif");
//selectWindow("Null_" + series+".tif");
//run("Split Channels");									
selectWindow(experiment+" - Null_" + series);
rename("Placeholder_Null.tif");
//selectWindow(experiment + " - Cell_" + series);
//rename("Placeholder_Cell.tif"); 

// Series

currentDirectory = getDirectory("Choose a Directory"); 
position_dir = currentDirectory + "Series" + series; 
File.makeDirectory(position_dir);

				
// RUN BLEACH CORRECTOIN 
run("Concatenate...", "  image1=Placeholder_Force.tif image2=Placeholder_Null.tif image3=[-- None --]");
run("Properties...", "channels=1 slices=1 frames=2");
rename("Bleach_Placeholder");

//run("Split Channels"); 
// Keep only the essential window 
//close("C2-Untitled"); 
Stack.setChannel(1);
//Stack.setFrame(1)

//run("Split Channels");
//selectWindow("C1-Bleach_Placeholder");
rename("Bleach_placehold");


run("Bleach Correction", "correction=[Histogram Matching]");
rename("Merged");

//close("C1-Bleach_Placeholder"); 
//selectWindow("DUP_C1-Bleach_Placeholder");

//run("Merge Channels...",  " c1=[DUP_C1-Bleach_Placeholder] c2=[C2-Bleach_Placeholder] create ignore keep") //c3=[C3-Bleach_Placeholder] c4=[C4-Bleach_Placeholder] create ignore keep"); //CHANGE 
//Stack.setDisplayMode("grayscale");


run("Make Substack...", "channels=1 frames=1"); //CHANGE with channel number
rename("Placeholder_Force.tif");

selectWindow("Merged")
run("Make Substack...", "channels=1 frames=2"); // CHANGE
rename("Null_Placeholder")

run("Concatenate...", "  image1=Null_Placeholder image2=Placeholder_Force.tif image3=[-- None --]");
close("\\Others")
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

//File.makeDirectory(position_dir + "/Trans");
//selectWindow("NullAppended_DeDrift");
//run("Make Substack...", "channels=2 frames=2");
//rename("Trans_Forced");

imagenum = 1;
saveAs(position_dir + "/Beads/Beads_Forced_" + IJ.pad(imagenum, 4) +".tif");
					
close("*")


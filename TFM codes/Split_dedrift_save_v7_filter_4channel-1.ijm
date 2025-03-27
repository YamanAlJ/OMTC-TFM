
run("Bio-Formats Macro Extensions");

userChosenDirectory = getDirectory("Choose a Directory");

processBioFormatFiles(userChosenDirectory);

function processBioFormatFiles(currentDirectory) {

	fileList = getFileList(currentDirectory);

	for (file = 0; file < fileList.length; file++) {
		Ext.isThisType(currentDirectory + fileList[file], supportedFileFormat);
		if (supportedFileFormat=="true") {
			Ext.setId(currentDirectory + fileList[file]);
			Ext.getSeriesCount(seriesCount);
			// Select your positions of interest for the analysis folder
			for (series = 9; series <= 13; series++) {
				//for (series = 1; series <= seriesCount; series++) {
				//record the Bio-Formats importer with the setup you need if different from below and change accordingly
				//   nOpen force and null image files then save as .tif
				run("Bio-Formats Importer", "open=[" + currentDirectory + fileList[file] + "] rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+series);
				rename("Placeholder_Force.tif");
				//saveAs("Tiff", "G:/My Drive/Master's/Research/Prostate/pyTFM/Cell Profiler Pipelines/ImageJ Pre-Processing/Placeholder_Force.tif");
				//  Update the Null image placements within the folder 
				run("Bio-Formats Importer", "open=[" + currentDirectory + fileList[file] + "] rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(series+29));
				rename("Placeholder_Null.tif");
				//run("Split Channels");
				//close("C2-Placeholder_Null.tif");
				//close("C3-Placeholder_Null.tif");
				//rename("Placeholder_Null.tif");

				//saveAs("Tiff", "G:/My Drive/Master's/Research/Prostate/pyTFM/Cell Profiler Pipelines/ImageJ Pre-Processing/Placeholder_Null.tif");
				// Add your macro code here
				// The following will create the necessary operations to dedrift and save the TIFF images
				// Step 1: Duplicate the null force as many times necessary to match that of the null force image. (Further manipulation may be required). 
				position_dir = currentDirectory + "Series" + series; 
				File.makeDirectory(position_dir);
				
				for (channum =1; channum <=4; channum++){ 
					selectWindow("Placeholder_Null.tif");
					run("Duplicate...", " title=" + channum);
				}

				//  Step 2: Merge the duplicate files together
				run("Merge Channels...", "c1=[1] c2=[2] c3=[3] c4=[4] create ignore keep");				
				// c3=[3] c4 = [4]"   list2[1] );
				rename("MergedNull");
				 // Step 3: Concatenate the files together for a global registration. Hyperstackreg uses the first frame as the reference anchor 
//				run("Concatenate...", "  image1=Placeholder_Force.tif image2=MergedNull image3=[-- None --]");
				// FOR PHOTOBLEACHING Beads
				run("Concatenate...", "  image1=Placeholder_Force.tif image2=MergedNull image3=[-- None --]");
				saveAs("TIFF", position_dir + "/Pre_bleach_placeholder.tif");
				rename("Bleach_Placeholder");
				Stack.setChannel(1);
				Stack.setFrame(1);
				run("Split Channels");
				selectWindow("C1-Bleach_Placeholder");
				//run("Enhance Contrast", "saturated=0.35");
				run("Bleach Correction", "correction=[Histogram Matching]");
				//run("Bandpass Filter...", "filter_large=30 filter_small=4 suppress=None tolerance=5 autoscale saturate process");
				close("C1-Bleach_Placeholder"); 
				selectWindow("DUP_C1-Bleach_Placeholder");
				run("Merge Channels...",  " c1=[DUP_C1-Bleach_Placeholder] c2=[C2-Bleach_Placeholder] c3=[C3-Bleach_Placeholder] c4=[C4-Bleach_Placeholder] create ignore keep"); //CHANGE 
				//saveAs("TIFF", position_dir + "/Post_bleach_placeholder.tif");
				SliceNum = nSlices;
				rename("Bleach_Corrected_Beads_Null_end"); 
				run("Make Substack...", "channels=1-4 frames=1-" + SliceNum/4-1); //CHANGE with channel number
				rename("Placeholder_Force.tif");
				selectWindow("Bleach_Corrected_Beads_Null_end");
				run("Make Substack...", "channels=1-4 frames=" + (SliceNum/4)); // CHANGE
				//saveAs("TIFF", position_dir + "/isthisreal.tif");
				rename("MergedNull"); 

				run("Concatenate...", "  image1=MergedNull image2=Placeholder_Force.tif image3=[-- None --]");
//				rename("NullAppended");
//				
//				// Split channels
//				title = getTitle(); 
//				run("Split Channels");
//				
//				selectWindow("C1-" + title);
//				run("StackReg ", "transformation=Translation");
//				Stack.setChannel(1);
//				Stack.setFrame(2);
//				run("Enhance Contrast", "saturated=0.35");
//				//run("Apply LUT");
//				run("Stack Contrast Adjustment", "is");
//				rename("NullAppended");
//				run("Merge Channels...", " c1=[C1-" + title +"] c2=[C2-" + title +"] c3=[C3-" + title +"] create ignore keep");	
//				
				rename("NullAppended");
				saveAs("TIFF", position_dir + "/NullAppended_uncorrected.tif");
				Stack.setDisplayMode("grayscale");
				//Step 4: Select correct window and channel reference for global registration 
				Stack.setChannel(1);
				Stack.setFrame(1);
//				run("HyperStackReg ", "transformation=Translation channel1 show");
//
				run("Correct 3D drift", "channel=1 multi_time_scale sub_pixel edge_enhance only=0 lowest=1 highest=1 max_shift_x=50 max_shift_y=50 max_shift_z=0");
				//run("Correct 3D drift", "channel=1 sub_pixel edge_enhance only=0 lowest=1 highest=1 max_shift_x=50 max_shift_y=50 max_shift_z=0");
				saveAs("TIFF", position_dir + "/DeDrift_uncropped.tif");
				makeRectangle(0,0, 1024, 1024);
				run("Crop");
				rename("NullAppended_DeDrift");
				close("\\Others");
				Stack.setDisplayMode("grayscale");
				// Crop to remove edge effects due to de-drift
				makeRectangle(50, 50,924, 924);
				run("Crop");
				//save();
				saveAs("TIFF", position_dir + "/NullAppended_DeDrift.tif");
				rename("NullAppended_DeDrift"); 
				
				// Step 5: Remove the null force image and reduce into a single image. 
				selectWindow("NullAppended_DeDrift"); 
				run("Make Substack...", "channels=1 frames=1");
				rename("Null Image");
				saveAs("Tiff", position_dir + "/Null.tif"); 
				close();
				// Step 6: Split tifs and resave for individual channels 
				// Beads 
				File.makeDirectory(position_dir + "/Beads");
				selectWindow("NullAppended_DeDrift"); 
				SliceNum = nSlices;
				print(SliceNum);
				run("Make Substack...", "channels=1 frames=2-" + SliceNum/4);
				rename("Beads_Forced");
				run("Stack to Images"); 
				//
				for (imagenum = 1; imagenum <= SliceNum/4-1; imagenum++){
					print("Beads_Forced-" + IJ.pad(imagenum, 4));
					selectWindow("Beads_Forced-" + IJ.pad(imagenum, 4));
					saveAs(position_dir + "/Beads/Beads_Forced_" + IJ.pad(imagenum, 4) +".tif");
					// close();
				} 
				selectWindow("NullAppended_DeDrift"); 
				close("\\Others");
				// Transmission  
				File.makeDirectory(position_dir + "/Trans");
				SliceNum = nSlices;
				print(SliceNum);
				run("Make Substack...", "channels=4 frames=2-" + SliceNum/4);
				rename("Trans_Forced");
				run("Stack to Images"); 
				for (imagenum = 1; imagenum <= (SliceNum/4-1); imagenum++){
					print("Trans_Forced-" + IJ.pad(imagenum, 4));
					selectWindow("Trans_Forced-" + IJ.pad(imagenum, 4));
					saveAs(position_dir + "/Trans/Trans_Forced_" + IJ.pad(imagenum, 4) +".tif");
				}

				// CellFluo
				selectWindow("NullAppended_DeDrift"); 
				close("\\Others");
				// Fluorescent Cells  
				File.makeDirectory(position_dir + "/CellFluo");
				SliceNum = nSlices;
				print(SliceNum);
				run("Make Substack...", "channels=3 frames=2-" + SliceNum/4);
				rename("CellFluo_Forced");
				run("Stack to Images"); 
				for (imagenum = 1; imagenum <= (SliceNum/4-1); imagenum++){
					print("CellFluo_Forced-" + IJ.pad(imagenum, 4));
					selectWindow("CellFluo_Forced-" + IJ.pad(imagenum, 4));
					saveAs(position_dir + "/CellFluo/CellFluo_Forced_" + IJ.pad(imagenum, 4) +".tif");
				}

	// Nucleus
				selectWindow("NullAppended_DeDrift"); 
				close("\\Others");
				// Nuclei 
				File.makeDirectory(position_dir + "/Nucleus");
				SliceNum = nSlices;
				print(SliceNum);
				run("Make Substack...", "channels=2 frames=2-" + SliceNum/4);
				rename("Nucleus_Forced");
				run("Stack to Images"); 
				for (imagenum = 1; imagenum <= (SliceNum/4-1); imagenum++){
					print("Nucleus_Forced-" + IJ.pad(imagenum, 4));
					selectWindow("Nucleus_Forced-" + IJ.pad(imagenum, 4));
					saveAs(position_dir + "/Nucleus/Nucleus_Forced_" + IJ.pad(imagenum, 4) +".tif");
				}

				selectWindow("NullAppended_DeDrift"); 
				close("*");

	

				// Step 6:  Separate each channel into a different file, name each accordingly. 
				//setOption("ExpandableArrays", true);
				//chanlist = newArray; 
				//chanlist[0] = "Beads" ;
				//chanlist[1] = "Fluocells";
				
				//chanlist[2] = "AR" 
				//chanlist[3] = "Nucleus" 
				// chanlist[4] = "
				//Array.print(chanlist)
				//Array.length()
				//for (i = 1, i <= chanlist.length, i++)  
				//Position_Dir = ; 
				//file.makeDirectory(splitDir); 
			}
			
		} else if (endsWith(fileList[file], "/")) {
			//processBioFormatFiles(currentDirectory + fileList[file]);
		}
	}
}

	//SliceNum = nSlices;
				//setSlice(nSlices); 
				//Stack.setChannel(1);
				//setSlice(1);
				//Stack.setChannel(1);

function saveroiarray(roi_array,directory_sereis){
	
	Array.print(roi_array);
	
	roiManager("Select", roi_array);
	
	roiManager("Save", directory_series +"/RoiSet.zip");
	
}

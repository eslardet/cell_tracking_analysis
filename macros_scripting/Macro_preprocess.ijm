

//folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/data/"
folder = "/Volumes/T7 Shield/Incucyte_data/"

//experiment = "107"
plate = "VID1723"
phase = "red"
pic = "1"

cells = newArray("D");	

// To remove any initial frame
start_frame = 1

column_start = 10
column_end = 10


for (i=0; i<cells.length; i++) {
	
	for (j=column_start; j<column_end+1; j++){
		
		well = cells[i] + d2s(j,0);
		
		path = folder + "raw_data/" + plate + "/" + plate + "_" + phase + "_" + well + "_" + pic + ".tif";
//		save_path = path;
		save_path = folder + "processed_data/" + plate + "/"  + plate + "_" + phase + "_" + well + "_" + pic + ".tif";
		
		open(path);
//		run("Slice Remover", "first=1 last=" + start_frame-1);
		
		for (n=1; n<=nSlices; n++) {
			setSlice(n);
			run("Enhance Contrast", "saturated=0.35");
			getMinAndMax(min, max);
			if (max > 10000) {
				setMinAndMax(10, 20);
			}
			run("Apply LUT", "slice");
		}
		save(save_path);
		close("*");
	}	
}



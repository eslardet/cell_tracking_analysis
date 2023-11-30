


//experiment = "109"
plate = "1713"

input = "/Volumes/T7 Shield/Incucyte_data/raw_data/" + plate + "/"
output = "/Volumes/T7 Shield/Incucyte_data/processed_data/" + plate + "/"

// To remove any initial frame
start_frame = 1

list = getFileList(input);
for (i=0; i<list.length; i++) {
	path = input + list[i];
	save_path = output + list[i];
	
	if (File.exists(save_path) == 0) {
		open(path);
		
		for (n=1; n<=nSlices; n++) {
			setSlice(n);
			run("Enhance Contrast", "saturated=0.35");
			getMinAndMax(min, max);
			if (max > 10000) {
				setMinAndMax(10, 20);
			}
			run("Apply LUT", "slice");
			getMinAndMax(min, max);
			setMinAndMax(60000, max);
			run("Apply LUT", "slice");
		}
		run("Slice Remover", "first=1 last=1 increment=1");
		save(save_path);
		
		close("*");
	}
}



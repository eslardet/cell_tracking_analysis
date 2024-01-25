

//experiment = "109"
plate = "1727"

input = "/Volumes/T7 Shield/Incucyte_data/processed_data/" + plate + "/"
output = "/Volumes/T7 Shield/Incucyte_data/processed_data_sliced/" + plate + "/"
if (File.exists(output) == 0) {
	File.makeDirectory(output);
}


// To remove any initial frame
end_frame = 216;

list = getFileList(input);
for (i=0; i<list.length; i++) {
	path = input + list[i];
	
	save_path = output + list[i];
	
	if (File.exists(save_path) == 0) {
		open(path);
		run("Slice Keeper", "first=2 last=" + end_frame + " increment=1");
		
		save(save_path);
		close("*");
	}
}



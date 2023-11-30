

//folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/data/"
folder = "/Volumes/T7 Shield/Incucyte_data/"

//experiment = "108"
plate = "VID1727"
pic = "1"

well_list = newArray("B3");	



for (i=0; i<well_list.length; i++) {
		
//	well = cells[i] + d2s(j,0);
	well = well_list[i];
	
	green_path = folder + "processed_data/" + plate + "/" + plate + "_green_" + well + "_" + pic + ".tif";
	red_path = folder + "processed_data/" + plate + "/" + plate + "_red_" + well + "_" + pic + ".tif";
//		save_path = path;
//		save_path = "/Users/el2021/Library/CloudStorage/OneDrive-ImperialCollegeLondon/PhD/Incucyte/animations/" + plate + "_" + well + ".avi";
	save_path = folder + "animations/" + plate + "_" + well + ".avi";
	
	save_ani = 1;
	if (File.exists(red_path) == 1) {
		open(red_path);
		image_red = getTitle();
		if (File.exists(green_path) == 1) {
			open(green_path);
			image_green = getTitle();
			run("Merge Channels...", "c1=["+image_red+ "] c2=["+image_green+ "] create");
		}
		else {
			run("Red");
		}
	}
	else {
		if (File.exists(green_path) == 1) {
			open(green_path);
			run("Green");
		}
		else {
			print("no green or red phase for well " + well);
			save_ani = 0;
		}
	}
	
	if (save_ani == 1) {
		setForegroundColor(255, 255, 0);
		run("Label...", "format=00:00:00 starting=0 interval=20 x=5 y=20 font=18 text=[] range=1-216 use");
		run("AVI... ", "compression=JPEG frame=10 save=save_path");
		close("*");
	}
}



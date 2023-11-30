

//folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/data/"
folder = "/Volumes/T7 Shield/Incucyte_data/"

//experiment = "108"
//plate = "VID1727"

pic = "1"

cells = newArray("C", "D");	
plate_list = newArray("1713");


// To remove any initial frame
start_frame = 1

column_start = 2
column_end = 11

for (p=0; p<plate_list.length; p++) {
	plate = plate_list[p];
	for (i=0; i<cells.length; i++) {
		
		for (j=column_start; j<column_end+1; j++){
			
			well = cells[i] + d2s(j,0);
			
			green_path = folder + "processed_data/" + plate + "/" + "VID" + plate + "_green_" + well + "_" + pic + ".tif";
			red_path = folder + "processed_data/" + plate + "/" + "VID" + plate + "_red_" + well + "_" + pic + ".tif";
			save_path = folder + "animations/" + plate + "/" + "VID" + plate + "_" + well + ".avi";
			
			save_ani = 1;
			if (File.exists(save_path) == 0) {
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
			}
			else {
				print("Already animated well " + well + "!");
				save_ani = 0;
			}
			
			if (save_ani == 1) {
				setForegroundColor(255, 255, 0);
				run("Label...", "format=00:00:00 starting=60 interval=60 x=5 y=20 font=18 text=[] range=1-24 use");
				run("AVI... ", "compression=JPEG frame=10 save=save_path");
				close("*");
			}
		}	
	}
}



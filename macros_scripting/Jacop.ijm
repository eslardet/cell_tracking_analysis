//run("TIFF Virtual Stack...", "open=[/Volumes/T7 Shield/Incucyte_data/processed_data/1737/VID1737_green_F4_1.tif]");
//run("Slice Keeper", "first=217 last=217 increment=1");
run("TIFF Virtual Stack...", "open=[/Volumes/T7 Shield/Incucyte_data/processed_data/1737/VID1737_red_F4_1.tif]");
run("Slice Keeper", "first=217 last=217 increment=1");
selectImage("VID1737_red_F4_1.tif");
run("Slice Keeper", "first=217 last=217 increment=1");
selectImage("VID1737_red_F4_1.tif");
close();
//selectImage("VID1737_green_F4_1.tif");
//close();
//run("JACoP ");
//
//selectImage("VID1737_red_F4_1.tif kept stack");
//selectImage("VID1737_red_F4_1.tif kept stack");
//run("JACoP ", "imga=[VID1737_green_F4_1.tif kept stack] imgb=[VID1737_red_F4_1.tif kept stack] thra=24152 thrb=10439 pearson overlap mm");

//dir=getDirectory("Select the directory");
//saveAs("Text", dir+"coloc");
//close("*");
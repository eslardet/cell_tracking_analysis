//inputFolder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/data"
//outputFolder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml"
//plate = "VID1711"
//phase = "red"
//pic = "1-1"
//
//t0 = time.time()
//for i in [10,11]:
//	for cell in ['B', 'C']:
//		well = cell + str(i)

folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/data/"

plate = "VID1713"
phase = "red"
//well = "B8"
pic = "1"
//cell = "B"


//numbers = newArray(1,2,3);
cells = newArray("B", "C", "D");	


for (i=0; i<cells.length; i++) {
	
//	for (j=0; j<numbers.length; j++){
	for (j=2; j<12; j++){
		
//		well = cells[i] + d2s(numbers[j],0);
		well = cells[i] + d2s(j,0);
		
		path = folder + plate + "_" + phase + "_" + well + "_" + pic + ".tif";
		path2 = folder + plate + "_" + phase + "_" + well + "_" + pic + "-1.tif";
		open(path);

		for (n=1; n<=nSlices; n++) {
			setSlice(n);
			run("Enhance Contrast", "saturated=0.35");
			run("Apply LUT", "slice");
		}
		
		save(path);
		close("*");	
	}	
}



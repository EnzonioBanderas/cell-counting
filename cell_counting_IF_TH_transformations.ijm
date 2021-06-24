dir = getDirectory("Choose a Transformation Directory"); //%%%//
print(dir);
files = getFileList(dir); //%%%//
print(files.length);
for(i=0;i<files.length;i++){
	//print(files[i]);
	if (endsWith (files[i], ".tif"))
	
	   {
		
		shortname = substring(files[i],0,indexOf(files[i],"."));
		print(shortname);
		main = dir+files[i];
		open(main);
	    
	   // run("Threshold...");
	    //setAutoThreshold("MaxEntropy dark");
		//setThreshold(161, 175);

		run("Split Channels");
//		files_split_string = split(files[i], ".");
		close(files[i]+ " (blue)");
		close(files[i]+ " (red)");
		main_green = files[i]+ " (green)";
        selectWindow(main_green);
//		run("Threshold...");
//		MaxEntropy RenyiEntropy Shanbhag Yen


		setAutoThreshold("MaxEntropy dark");
//		setAutoThreshold("RenyiEntropy dark");
//		setAutoThreshold("Shanbhag dark");
//		setAutoThreshold("Yen dark");

//		waitForUser("set threshold");
		run("Convert to Mask");
		TH_masks_dir = dir+"TH_masks";
		File.makeDirectory(TH_masks_dir);
		saveAs("tiff", TH_masks_dir+File.separator+shortname+"_TH_mask.tif");
		
//		run("Set Measurements...", "centroid nan redirect=None decimal=0");
//        run("Analyze Particles...", "size=0-50 pixel circularity=0-1.00 Include holes display add");
//	        
//		saveAs("Results", dir+ shortname+ ".csv");
//		run("Clear Results");
//		roiManager("Delete");
//		print("\\Clear");

		}	

	run("Close All");
	
}
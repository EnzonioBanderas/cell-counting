dir = getDirectory("Choose a Directory");
files = getFileList(dir);

for(i=0;i<files.length;i++){
	shortname = substring(files[i],0,indexOf(files[i],"."));
	if (endsWith (files[i], ".tif"))
	
	   {
		
		open(dir+files[i]);
	    
	   // run("Threshold...");
	    //setAutoThreshold("MaxEntropy dark");
		//setThreshold(161, 175);

		run("Split Channels");
//		files_split_string = split(files[i], ".");
        selectWindow(files[i]+ " (red)");
//		run("Threshold...");
//		MaxEntropy RenyiEntropy Shanbhag Yen

		setAutoThreshold("MaxEntropy dark");
//		setAutoThreshold("RenyiEntropy dark");
//		setAutoThreshold("Shanbhag dark");
//		setAutoThreshold("Yen dark");

//		waitForUser("set threshold");
		
		run("Set Measurements...", "centroid nan redirect=None decimal=0");
        run("Analyze Particles...", "size=0-400 pixel circularity=0-1.00 Include holes display clear add");
		saveAs("Results", dir + shortname+ ".csv");


		}	

	run("Close All");
	
}


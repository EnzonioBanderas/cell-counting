dir = getDirectory("Choose Data Directory");
print(dir);
dirs = getFileList(dir); /////////////

//for(u=0;u<1.1;u++){	/////////////
for(u=0;u<dirs.length;u++){	/////////////
print(dirs[u]);
if (endsWith (dirs[u], File.separator)){ //////////////
// files = getFileList(dir); //%%%//
files = getFileList(dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator);
print(files.length);
//for(i=0;i<2.1;i++){
for(i=0;i<files.length;i++){
	//print(files[i]);
	if (endsWith (files[i], ".tif"))
	
	   {
		
		shortname = substring(files[i],0,indexOf(files[i],"."));
		print(shortname);
		structureStack_name = shortname+"_structureStack.tif";
		structureStack_path = dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator+"structure_stacks"+File.separator+structureStack_name;
		open(structureStack_path);
		main = dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator+files[i];
		open(main);
	    
	   // run("Threshold...");
	    //setAutoThreshold("MaxEntropy dark");
		//setThreshold(161, 175);

		run("Split Channels");
//		files_split_string = split(files[i], ".");
		close(files[i]+ " (blue)");
		close(files[i]+ " (green)");
		main_red = files[i]+ " (red)";
        selectWindow(main_red);
//		run("Threshold...");
//		MaxEntropy RenyiEntropy Shanbhag Yen

		selectWindow(structureStack_name);
		for (s = 1; s <= nSlices; s++) {
			selectWindow(structureStack_name);
		    setSlice(s);

			imageCalculator("AND create", main_red, structureStack_name);

	//		setAutoThreshold("MaxEntropy dark");
			run("Auto Threshold", "method=MaxEntropy ignore_black white");
	//		setAutoThreshold("RenyiEntropy dark");
	//		setAutoThreshold("Shanbhag dark");
	//		setAutoThreshold("Yen dark");
	
	//		waitForUser("set threshold");
			
			run("Set Measurements...", "centroid nan redirect=None decimal=0");
	        run("Analyze Particles...", "size=0-50 pixel circularity=0-1.00 Include holes display add");

	        close();
	        close("Exception");
		}
		saveAs("Results", dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator+ shortname+ ".csv");
		run("Clear Results");
		roiManager("Delete");
		print("\\Clear");

		}	

	run("Close All");
	
}
} ////////////////
} ////////////////
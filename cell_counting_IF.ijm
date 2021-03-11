// dir = getDirectory("Choose a Transformation Directory"); //%%%//
dir = getDirectory("Choose Data Directory");
print(dir);
dirs = getFileList(dir); /////////////

for(u=0;u<dirs.length;u++){	/////////////
print(dirs[u]);
if (endsWith (dirs[u], File.separator)){ //////////////
// files = getFileList(dir); //%%%//
files = getFileList(dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator);
print(files.length);
for(i=0;i<files.length;i++){
	//print(files[i]);
	if (endsWith (files[i], ".tif"))
	
	   {
		
		shortname = substring(files[i],0,indexOf(files[i],"."));
		print(shortname);
		open(dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator+files[i]);
	    
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
        run("Analyze Particles...", "size=0-50 pixel circularity=0-1.00 Include holes display clear add");
		saveAs("Results", dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator+ shortname+ ".csv");


		}	

	run("Close All");
	
}
} ////////////////
} ////////////////
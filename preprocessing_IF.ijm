dir = getDirectory("Choose a Directory");
files = getFileList(dir);

for(i=0;i<files.length;i++){
	shortname = substring(files[i],0,indexOf(files[i],"."));
	if (endsWith (files[i], "AF594.tif"))
	
	   {
		
		open(dir+files[i]);
//		run("Zoom", "Original Scale")
//		run("Channels Tool...");
//        run("Make Composite", "Display Mode=Composite");
//		run("Split Channels");
		selectWindow(files[i]);
		run("Rotate 90 Degrees Right");
		run("Rotate... "); //use the dialog window to adjust the angle if the slice is a bit tilted. IMP! click on Preview to see the changes. Press OK when done.
        setTool("Rectangle");
        waitForUser("Select region for cropping. \nClick ‘OK’ when done.");
        run("Crop");
        selectWindow(files[i]);
        run("Scale...", "x=- y=- width=1140 height=800 interpolation=Bilinear fill average create"); //to match reference atlas
        run("Coordinates...", "left=0 right=1140 top=0 bottom=800");
		saveAs(dir + "downsampled_" + shortname + ".tif");

		}	

	run("Close All");
	
}


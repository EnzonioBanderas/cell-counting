run("Record...")

dir = getDirectory("Choose a Original Directory with ParamTable.csv inside, exact match to ParamTable.csv");
files = getFileList(dir);

iTif = -1;
for(i=0;i<files.length;i++){
	if (endsWith (files[i], ".tif"))

	   {

		Table.open(dir+"ParamTable.csv");

		iTif = iTif + 1;
		ro_array = Table.getColumn("ro");
		ro = ro_array[iTif]; 
		mr1_array = Table.getColumn("mr1");
		mr1 = mr1_array[iTif]; 
		mr2_array = Table.getColumn("mr2");
		mr2 = mr2_array[iTif]; 
		mr3_array = Table.getColumn("mr3");
		mr3 = mr3_array[iTif]; 
		mr4_array = Table.getColumn("mr4");
		mr4 = mr4_array[iTif]; 

		print(ro);
		
		shortname = substring(files[i],0,indexOf(files[i],"."));
		open(dir+files[i]);
//		run("Zoom", "Original Scale")
//		run("Channels Tool...");
//        run("Make Composite", "Display Mode=Composite");
//		run("Split Channels");
		selectWindow(files[i]);
		run("Rotate 90 Degrees Right");
		run("View 100%");
		run("Original Scale");
//		run("Rotate... "); //use the dialog window to adjust the angle if the slice is a bit tilted. IMP! click on Preview to see the changes. Press OK when done.
        //waitForUser("Rotate. \nClick ‘OK’ when done.");
        run("Rotate... ", "angle="+ro+" grid=1 interpolation=Bilinear");
        setTool("Rectangle");
        makeRectangle(mr1, mr2, mr3, mr4);
        //waitForUser("Select region for cropping. \nClick ‘OK’ when done.");
        run("Crop");
        selectWindow(files[i]);
        run("Scale...", "x=- y=- width=1140 height=800 interpolation=Bilinear fill average create"); //to match reference atlas
        run("Coordinates...", "left=0 right=1140 top=0 bottom=800");
		saveAs(dir + "downsampled_" + shortname + ".tif");
		print(files[i]);

		}	

	run("Close All");
	
}
//selectWindow("Recorder");
//save(dir+"Recorder.txt");
selectWindow("Log");
save(dir+"Log.txt");


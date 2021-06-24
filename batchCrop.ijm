// Batchcrop
directory = "/home/enzo/Pictures/"
filelist = getFileList(directory) ;
for (i = 0; i < lengthOf(filelist); i++) {
    if (endsWith(filelist[i], ".png")) { 
    	print(filelist[i]);
        //run("Open", directory + File.separator + filelist[i]);
        open(filelist[i]);
        makeRectangle(2389, 190, 4071-2389, 1374-198);
        //run("Crop");
        run("Crop");
        save(filelist[i]+'_cropped.png');
    } 
}
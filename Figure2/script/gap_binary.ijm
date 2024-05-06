close("*"); close("Results"); close("ROI Manager");

// set the folder containing the raw images to be processed
folder = "ftw_paper_figures_v2/fig2_gap";
Table.open("D:/"+folder+"/data/parameters.csv");

//for (i = 0; i < Table.size; i++) { //Table.size
	i = 4;
	// set parameters' value
	well = Table.get("wells", i);
	frame_strt_global = Table.get("frame_strt_global", i);
	frame_scratch = Table.get("frame_scratch", i);
	frame_stop_global = Table.get("frame_stop_global", i);
	angle = Table.get("angle",i);
	x = Table.get("x_local",i);
	y = Table.get("y_local",i);
	frame_strt_local = Table.get("frame_strt_local", i);
	frame_stop_local = Table.get("frame_stop_local", i);
	pxSize = Table.get("pxSize",i);
		
	// load images
	open("D:/"+folder+"/img/s"+well+"_ros_ff1_10_aligned_cropped.tif"); rename("ros");
	
	// deal with the black frames in the well 89 and 30
	if (well==89) {
		selectWindow("ros");setSlice(4);run("Delete Slice");
		selectWindow("ros");setSlice(7);run("Delete Slice");
	}
	if (well==30) {
		selectWindow("ros");setSlice(6);run("Delete Slice");
		selectWindow("ros");setSlice(9);run("Delete Slice");
	}
	
	// get the third frame of local view
	run("Duplicate...", "duplicate range="+frame_strt_global+"-"+frame_stop_global);
	run("Duplicate...", "duplicate range="+frame_strt_local+"-"+frame_stop_local);
	setSlice(3);
	run("Duplicate...", " ");	
	run("Rotate... ", "angle="+angle+" grid=1 interpolation=Bilinear stack");
	makeRectangle(x, y, 201, Math.round(2000/pxSize));	
	run("Crop");
	
	// save image
	saveAs("Tiff", "D:/"+folder+"/img/gap/s"+well+"_ros_gap.tif");
		
	close("*");	
}

/*
close("*"); close("Results"); close("ROI Manager");

// set the folder containing the raw images to be processed
folder = "ftw_paper_figures_v2/fig2_gap";
wells = newArray(48,68,55,87,89,12,18,30);
signals = newArray("cy5","ros","dic");

// angle of rotating images
angle = newArray(-129.584,-133.108,-124.173,-141.576,-129.588,-135.51,-138.421,-136.831);

// coordinate of upper left corner for local view
xc2 = newArray(1008,867,1080,1548,1350,1548,1648,1200);
yc2 = newArray( 356,186, 824,  98, 776, 869, 502, 434);

// start and stop frames of movies
frame_strt_global = newArray(4,4,4,4,6,6,6,5);
frame_stop_global = newArray(20,20,21,20,15,22,21,22);

// frame of scratch, adjust the display range for that frame
frame_scratch = newArray(8,8,8,8,7,11,11,9);

// start and stop frames of local view
frame_strt_local = newArray( 5, 6, 5, 7,3, 7, 8, 7);
frame_stop_local = newArray(11,12,11,13,9,13,14,13);

// for each well, get the third frame for gap size measurement
for (i = 0; i < wells.length; i++) {
	
	open("D:/"+folder+"/img/s"+wells[i]+"_"+signals[1]+"_ff1_10_aligned_cropped.tif"); rename("ros");
	
	if (wells[i]==89) {
		selectWindow("ros");setSlice(4);run("Delete Slice");
		selectWindow("ros");setSlice(7);run("Delete Slice");
	}
	if (wells[i]==30) {
		selectWindow("ros");setSlice(6);run("Delete Slice");
		selectWindow("ros");setSlice(9);run("Delete Slice");
	}
	
	// get the third frame of local view
	run("Duplicate...", "duplicate range="+frame_strt_global[i]+"-"+frame_stop_global[i]);
	run("Duplicate...", "duplicate range="+frame_strt_local[i]+"-"+frame_stop_local[i]);
	
	setSlice(3);
	run("Duplicate...", " ");
	
	run("Rotate... ", "angle="+angle[i]+" grid=1 interpolation=Bilinear stack");
	makeRectangle(xc2[i], yc2[i], 201, 1600);	
	run("Crop");
	saveAs("Tiff", "D:/"+folder+"/img/s"+wells[i]+"_ros_gap.tif");
		
	close("*");
}
*/
close("*"); close("Results"); close("ROI Manager");

// set the folder containing the raw images to be processed
folder = "ftw_paper_figures_v2/fig2_gap";
Table.open("D:/"+folder+"/data/parameters.csv");

//for (i = 0; i < Table.size; i++) { //Table.size
	i = 0;
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
	open("D:/"+folder+"/img/s"+well+"_cy5_ff1_10_aligned_cropped.tif"); rename("cy5");
	open("D:/"+folder+"/img/s"+well+"_ros_ff1_10_aligned_cropped.tif"); rename("ros");
	
	// deal with the black frames in the well 89 and 30
	if (well==89) {
		selectWindow("cy5");setSlice(4);run("Delete Slice");
		selectWindow("cy5");setSlice(7);run("Delete Slice");
		selectWindow("ros");setSlice(4);run("Delete Slice");
		selectWindow("ros");setSlice(7);run("Delete Slice");
	}
	if (well==30) {
		selectWindow("cy5");setSlice(6);run("Delete Slice");
		selectWindow("cy5");setSlice(9);run("Delete Slice");
		selectWindow("ros");setSlice(6);run("Delete Slice");
		selectWindow("ros");setSlice(9);run("Delete Slice");
	}
	
	// get differential ROS
	selectWindow("ros");
	run("Duplicate...", "duplicate range=1-"+(nSlices-1));
	selectWindow("ros");
	run("Duplicate...", " ");
	run("Concatenate...", "open image1=ros-2 image2=ros-1 image3=[-- None --]");
	imageCalculator("Subtract create stack", "ros","Untitled");
	rename("dros");
	
	// composite images
	run("Merge Channels...", "c1=cy5 c2=dros create");
	run("Median...", "radius=2");
	setSlice(1); run("Cyan Hot"); setMinAndMax(350, 600);
	setSlice(2); run("Yellow"); setMinAndMax(30, 90);
	
	selectWindow("Merged");
	run("Duplicate...", "duplicate frames="+frame_strt_global+"-"+(frame_scratch-1));
	run("RGB Color", "frames keep");
	rename("before");

	selectWindow("Merged");
	run("Duplicate...", "duplicate frames="+(frame_scratch+1)+"-"+frame_stop_global);
	run("RGB Color", "frames keep");
	rename("after");

	selectWindow("Merged");
	run("Duplicate...", "duplicate frames="+frame_scratch);
	setSlice(2);setMinAndMax(60, 100);
	run("RGB Color", "frames keep");
	rename("scratch");
	
	// export the movie of global view
	run("Concatenate...", "open image1=before image2=scratch image3=after image4=[-- None --]");
	run("AVI... ", "compression=JPEG frame=7 save=D:/"+folder+"/fig/s"+well+"_global.avi");
	
	// get local view and export the montage of them and the movie
	run("Rotate... ", "angle="+angle+" grid=1 interpolation=Bilinear stack");
	makeRectangle(x, y, 201, Math.round(2000/pxSize)); // 	
	run("Crop");
	run("Duplicate...", "duplicate range="+frame_strt_local+"-"+frame_stop_local);
	run("Make Montage...", "columns=7 rows=1 scale=1 border=6");
	run("Rotate 90 Degrees Right");
	saveAs("Tiff", "D:/"+folder+"/fig/s"+well+"_composite_med2_montage7_wo_scaleBar.tif");
	
	if (i==0) {		
		run("Scale Bar...", "width=400 height=24 font=100 color=White background=None location=[Upper Right] bold");
		saveAs("Tiff", "D:/"+folder+"/fig/s"+well+"_composite_med2_montage7_w_scaleBar.tif");
	}
	
	selectWindow("Untitled-1");
	run("AVI... ", "compression=JPEG frame=7 save=D:/"+folder+"/fig/s"+well+"_local.avi");
	
	close("*");
	
}

/*
 * makeRectangle(1008, 356, 201, 1580);

wells  = newArray(  48,  68,  55,87,89,12,18,30);
pxSize = newArray(1.26,1.26,1.26);
signals = newArray("cy5","ros","dic");

// angle of rotating images
angle = newArray(-129.584,-133.108,-124.173,-141.576,-129.588,-135.51,-138.421,-136.831);

// coordinate of upper left corner for local view
xc2 = newArray(1008,867,1080,1548,1350,1548,1648,1200);
yc2 = newArray( 356,186, 824,  98, 776, 869, 502, 434);

// start and stop frames of global view
frame_strt_global = newArray(4,4,4,4,6,6,6,5);
frame_stop_global = newArray(20,20,21,20,15,22,21,22);

// frame of scratch, adjust the display range for that frame
frame_scratch = newArray(8,8,8,8,7,11,11,9);

// start and stop frames of local view
frame_strt_local = newArray( 5, 6, 5, 7,3, 7, 8, 7);
frame_stop_local = newArray(11,12,11,13,9,13,14,13);

// for each well, generate 2 movies (global and local) and a montage of local view
for (i = 0; i < wells.length; i++) {
	
	// load images
	open("D:/"+folder+"/img/s"+wells[i]+"_"+signals[0]+"_ff1_10_aligned_cropped.tif"); rename("cy5");
	open("D:/"+folder+"/img/s"+wells[i]+"_"+signals[1]+"_ff1_10_aligned_cropped.tif"); rename("ros");
	
	// deal with the black frames in the well 89 and 30
	if (wells[i]==89) {
		selectWindow("cy5");setSlice(4);run("Delete Slice");
		selectWindow("cy5");setSlice(7);run("Delete Slice");
		selectWindow("ros");setSlice(4);run("Delete Slice");
		selectWindow("ros");setSlice(7);run("Delete Slice");
	}
	if (wells[i]==30) {
		selectWindow("cy5");setSlice(6);run("Delete Slice");
		selectWindow("cy5");setSlice(9);run("Delete Slice");
		selectWindow("ros");setSlice(6);run("Delete Slice");
		selectWindow("ros");setSlice(9);run("Delete Slice");
	}
	
	// get differential ROS
	selectWindow("ros");
	run("Duplicate...", "duplicate range=1-"+(nSlices-1));
	selectWindow("ros");
	run("Duplicate...", " ");
	run("Concatenate...", "open image1=ros-2 image2=ros-1 image3=[-- None --]");
	imageCalculator("Subtract create stack", "ros","Untitled");
	rename("dros");
	
	// composite images
	run("Merge Channels...", "c1=cy5 c2=dros create");
	run("Median...", "radius=2");
	setSlice(1); run("Cyan Hot"); setMinAndMax(350, 600);
	setSlice(2); run("Yellow"); setMinAndMax(30, 90);
	
	selectWindow("Merged");
	run("Duplicate...", "duplicate frames="+frame_strt_global[i]+"-"+(frame_scratch[i]-1));
	run("RGB Color", "frames keep");
	rename("before");

	selectWindow("Merged");
	run("Duplicate...", "duplicate frames="+(frame_scratch[i]+1)+"-"+frame_stop_global[i]);
	run("RGB Color", "frames keep");
	rename("after");

	selectWindow("Merged");
	run("Duplicate...", "duplicate frames="+frame_scratch[i]);
	setSlice(2);setMinAndMax(60, 100);
	run("RGB Color", "frames keep");
	rename("scratch");
	
	// export the movie of global view
	run("Concatenate...", "open image1=before image2=scratch image3=after image4=[-- None --]");
	run("AVI... ", "compression=JPEG frame=7 save=D:/"+folder+"/fig/s"+wells[i]+"_global.avi");
	
	// get local view and export the montage of them and the movie
	run("Rotate... ", "angle="+angle[i]+" grid=1 interpolation=Bilinear stack");
	makeRectangle(xc2[i], yc2[i], 201, 1600);
	run("Crop");
	run("Duplicate...", "duplicate range="+frame_strt_local[i]+"-"+frame_stop_local[i]);
	run("Make Montage...", "columns=7 rows=1 scale=1 border=6");
	run("Rotate 90 Degrees Right");
	saveAs("Tiff", "D:/"+folder+"/fig/s"+wells[i]+"_composite_med2_montage7_wo_scaleBar.tif");
	
	if (i==0) {		
		run("Scale Bar...", "width=400 height=24 font=100 color=White background=None location=[Upper Right] bold");
		saveAs("Tiff", "D:/"+folder+"/fig/s"+wells[i]+"_composite_med2_montage7_w_scaleBar.tif");
	}
	
	selectWindow("Untitled-1");
	run("AVI... ", "compression=JPEG frame=7 save=D:/"+folder+"/fig/s"+wells[i]+"_local.avi");
	
	close("*");
}

*/
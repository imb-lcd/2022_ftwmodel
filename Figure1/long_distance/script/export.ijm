close("*"); close("Results"); close("ROI Manager");

// set the folder containing the images to be exported
folder = "ftw_paper_figures_v2/fig1_long_distance";

// load the image to be exported
open("D:/"+folder+"/img/s19_cy5_ff1_10.tif");
rename("cy5");

// determine display range
makeRectangle(992, 2527, 45, 45);
run("Measure");
meanLv = getResult("Mean", 0);
print(meanLv);
close("Results");
run("Select None");
setMinAndMax(meanLv, meanLv+255);
run("Cyan Hot");

// apply median filter to reduce noise
run("Median...", "radius=2 stack");

// convert to RGB images and save
setSlice(13);
run("Duplicate...", " ");
makeRectangle(0, 500, 4500, 4500);
run("Crop");
run("RGB Color");
//saveAs("Tiff", "D:/"+folder+"/fig/s19_cy5_ff1_10_med2_crop4500_wo_scaleBar.tif");
// add a scale bar
run("Scale Bar...", "width=1000 height=36 font=150 color=White background=None location=[Upper Right] bold");
//saveAs("Tiff", "D:/"+folder+"/fig/s19_cy5_ff1_10_med2_crop4500_w_scaleBar.tif");
close("s19*");

// export the images around the initiation
selectWindow("cy5");
run("Duplicate...", "duplicate range=1-5");
//makeRectangle(212, 3813, 901, 901);
makeRectangle(653-450, 4294-450, 901, 901);
makeRectangle(199, 3845, 901, 901);
run("Crop");
run("RGB Color");
run("Image Sequence... ", "select=D:/"+folder+"/fig/ dir=D:/"+folder+"/fig/ format=TIFF name=s19_cy5_ff1_10_crop900_f digits=1");

setSlice(1);
run("Duplicate...", " ");
run("Scale Bar...", "width=400 height=12 font=50 color=White background=None location=[Upper Right] bold");
saveAs("Tiff", "D:/"+folder+"/fig/s19_cy5_ff1_10_crop900_f1_w_scaleBar.tif");

// export the bright field images around the initiation
open("D:/"+folder+"/img/s19_dic_ff1_10.tif");
rename("bright_field")
//makeRectangle(212, 3813, 901, 901);
makeRectangle(653-450, 4294-450, 901, 901);
makeRectangle(199, 3845, 901, 901);
run("Crop");
run("Enhance Contrast", "saturated=2");
run("Duplicate...", "duplicate range=1-5");
run("RGB Color");
run("Image Sequence... ", "select=D:/"+folder+"/fig/ dir=D:/"+folder+"/fig/ format=TIFF name=s19_dic_ff1_10_crop900_f digits=1");

// export cropped images
folder = "ftw_paper_figures_v2/fig1_long_distance";
open("D:/"+folder+"/img/s19_cy5_ff1_10_sliced.tif");
run("Median...", "radius=2 stack");
setMinAndMax(234, 526);
run("Cyan Hot");
run("Make Montage...", "columns=16 rows=1 scale=1 first=1 last=16 border=12");
run("RGB Color");
saveAs("Tiff", "D:/"+folder+"/fig/s19_cy5_ff1_10_sliced_montaged.tif");

// movie
selectWindow("cy5");
run("Duplicate...", "duplicate range=1-27");
makeRectangle(0, 500, 4500, 4500);
run("Crop");
run("RGB Color");

setSlice(1);
run("Duplicate...", " ");
setForegroundColor(0, 255, 255);
makeOval(465, 3582, 394, 394);
run("Dotted Line", "line=20 dash=20,20");

run("Concatenate...", "open image1=cy5-3 image2=cy5-2 image3=[-- None --]");
setSlice(2);
run("Delete Slice");

run("Scale Bar...", "width=1000 height=36 font=150 color=White background=None location=[Upper Right] bold label");
run("AVI... ", "compression=JPEG frame=3 save=D:/ftw_paper_figures_v2/fig1_long_distance/fig/s19_cy5_ff1_10_med2_crop4500.avi");

// try to use mm for scale bar
//Stack.setXUnit("mm");
//Stack.setYUnit("mm");
//run("Properties...", "channels=1 slices=1 frames=27 pixel_width=0.0012660 pixel_height=0.0012660 voxel_depth=0.0013000");
// need manually add scale bar (why?)
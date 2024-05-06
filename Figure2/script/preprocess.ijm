close("*"); close("Results"); close("ROI Manager");

// set the folder containing the raw images to be processed
folder = "ftw_paper_figures_v2/fig2_gap";
Table.open("D:/"+folder+"/data/preprocess.csv");
wells   = newArray(48,48,48,68,68,68,55,55,55,87,87,87,89,89,89,12,12,12,18,18,18,30,30,30);
ch      = newArray(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3);
signals = newArray("cy5","ros","dic","cy5","ros","dic","cy5","ros","dic","cy5","ros","dic","cy5","ros","dic","ros","cy5","dic","ros","cy5","dic","ros","cy5","dic");

// coordinate of upper left corner for cropping aligned images
xc = newArray(132,132,132,56,56,56,64,64,64,56,56,56,112,112,112,112,112,112,112,112,112,112,112,112);
yc = newArray(132,132,132,164,164,164,164,164,164,156,156,156,112,112,112,112,112,112,112,112,112,112,112,112);

// for each channel in each well, do flat-field correction and alignment
for (i = 0; i < wells.length; i++) {
	// set parameters
	fname = Table.getString("imgname", i); // filename of the image stack to be processed
	well = Table.get("well", i); // the well where the image stack is obtained
	ch = Table.get("ch", i); // the channel of the image stack
	signal = Table.getString("signal", i); // the signal of the channel	
	pxSize = Table.get("pxSize", i); // pixel size of the image stack
	xc = Table.get("xc",i);
	yc = Table.get("yc",i);
	
	// load raw image
	open("D:/"+folder+"/raw/s"+well+"c"+ch+"_ORG.tif");
	rename("source");
	
	// measure the mean intensity of the first frame
	setSlice(1); 
	run("Select All"); 
	run("Measure");
	meanLv = getResult("Mean", 0);	print(meanLv);	close("Results");	
	
	// background estimation using the first frame
	run("Duplicate...", " "); rename("background");
	run("Gaussian Blur...", "sigma=10");
	
	// flat-field correction based on the background from the first frame
	selectWindow("source");
	run("Calculator Plus", "i1=source i2=background "
	+"operation=[Divide: i2 = (i1/i2) x k1 + k2] k1="+meanLv+" k2=0 create");
	
	// set the properties of images
	Stack.setXUnit("um");
	Stack.setYUnit("um");
	run("Properties...", "channels=1 slices=1 frames="+nSlices
	+" pixel_width="+pxSize+" pixel_height="+pxSize+" voxel_depth=1.266");
	
	// alignment
	print("s"+well);
	run("MultiStackReg", "stack_1=Result action_1=[Load Transformation File]"+
	" file_1=[D:/"+folder+"/TransformationMatrices/TransformationMatrices_s"+well+".txt]"+
	" stack_2=None action_2=Ignore file_2=[] transformation=[Rigid Body]");
	makeRectangle(xc, yc, 2800, 2800);
	run("Crop");
	
	// save the processed images in the img folder
	saveAs("Tiff", "D:/"+folder+"/img/s"+well+"_"+signal+"_ff1_10_aligned_cropped.tif");
	
	close("*");

}

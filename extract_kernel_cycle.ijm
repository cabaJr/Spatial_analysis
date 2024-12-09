// Macro to define grid ROIs, save them to the ROI Manager, and use Multi Measure

// Define parameters
kernelSize = 5;  // Size of each kernel
pixelSize = 1;
inputDir = getDirectory("Select the input folder containing the images");
fileList = getFileList(inputDir);

// Loop through the file list to collect unique base filenames
for (i = 0; i < lengthOf(fileList); i++) {
    fileName = fileList[i];
msg1 = "extracting values from: "+ fileName;
print(msg1);
    // Only process .tif files
    if (endsWith(fileName, ".tif") | endsWith(fileName, ".tiff")) {
    	setBatchMode(true);
    	filepath = inputDir + File.separator + fileName;
    	open(filepath);
title = File.nameWithoutExtension();

setImgProperties(pixelSize);

width = 0;
height = 0;
channels = 0;
slices = 0;
frames = 0;
getDimensions(width, height, channels, slices, frames); 

//get directory where to save files
output = inputDir + title + "_results";

// create grid
createGrid(width, height, kernelSize);

// get centroids coordinates
getROIsCentroids(output);

// get intensity values
getROIsMean(output, 2, title);
getROIsMean(output, 1, title);

close(fileName);
setBatchMode(false);
msg = "Values extracted from: "+ title;
print(msg);
    } //end if statement
    } // end for cycle

function createGrid(imageWidth, imageHeight, kernelSize) { 
// function to create a grid of ROIs given a specified kernel size
	// Calculate the number of kernels per row and column
	kernelsPerRow = imageWidth / kernelSize;
	kernelsPerColumn = imageHeight / kernelSize;
	// Create grid of ROIs and add each to the ROI Manager
	roiManager("Reset");
	for (y = 0; y < kernelsPerColumn; y++) {
	    for (x = 0; x < kernelsPerRow; x++) {
	        // Calculate the top-left position for each kernel
	        xPos = x * kernelSize;
	        yPos = y * kernelSize;
	        // Define the ROI for the kernel and add it to the ROI Manager
	        makeRectangle(xPos, yPos, kernelSize, kernelSize);
	        roiManager("Add");
	    }
	}
}

function getROIsCentroids(dirpath) {
	// Get coordinates of each ROI centroid as X and Y values and save in csv file
	run("Set Measurements...", "centroid redirect=None decimal=3");
	roiManager("Deselect");
	roiManager("multi-measure");
	pathcsv = dirpath+"/"+"grid_centroids.csv";
	saveAs("Results", pathcsv);
	close("Results");
	}
	
function getROIsMean(dirpath, channel, title) {
	title = "Ch"+channel;
	run("Select All");
	run("Duplicate...", "title="+title+" duplicate channels="+channel);
	// extract intensity value from all ROIs from all frames
	run("Set Measurements...", "mean redirect=None decimal=3");
	roiManager("Multi Measure");
	pathcsv = dirpath+"/"+title+"_"+channel+"_grid_vals.csv";
	saveAs("Results", pathcsv);
	close("Results");
	close(title);
}

function setImgProperties(pxSize){
//	args
	// pxSize	pixelSize 
	
	Stack.setXUnit("um");
	run("Properties...", "pixel_width="+pxSize+" pixel_height="+pxSize+" voxel_depth="+pxSize+" frame=[1800 sec]");
	}
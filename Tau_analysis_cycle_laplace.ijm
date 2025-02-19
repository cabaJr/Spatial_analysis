/* Add workflow description

*/

// general parameters
pixelSize = 2.82;
run("Set Measurements...", "area mean standard min centroid center perimeter bounding shape integrated median stack limit redirect=None decimal=3");

// get channels of interest from user
channel_fred = getString("What is the channel containing FRED?", "1");	
channel_gcamp = getString("What is the channel containing GCaMP?", "2");

// start of do cycle
inputDir = getDirectory("Select the input folder containing the images");
fileList = getFileList(inputDir);
// Loop through the file list to collect unique base filenames
for (i = 0; i < lengthOf(fileList); i++) {
    fileName = fileList[i];

    // Only process .tif files
    if (endsWith(fileName, ".tif") | endsWith(fileName, ".tiff")) {
// open image, get ID, get title, set img properties, and create result folder
open(inputDir + fileName);
originalImg = getImageID();
title = File.nameWithoutExtension();
setImgProperties(pixelSize);

FileDir = File.directory;

// create result folder
resultDir = FileDir + title + "_results";
File.makeDirectory(resultDir);	

// create global variable for thresholded image
var thresholded_img = -1;
prefixwhole = "whole";
prefixscn = "SCN";

// get whole sliceROI from function
print("Defining slice area");
wholeSliceROI = defineSlice(channel_gcamp, 280, title, resultDir, originalImg);
print("Area defined succesfully!");

// threshold image
area = true;
gammasave = true;
closeImg = false;
threshold = 300;
filterSize = 4;
print("Preparing stack to be thresholded");
//threshold_path = threshold_prep(originalImg, channel_fred, area, wholeSliceROI, resultDir, gammasave, closeImg);
threshold_path = laplaceThresholdStack(originalImg, channel_fred, area, wholeSliceROI, resultDir, closeImg, threshold); 
print("Stack thresholding successful!");

// perform aggregated particle analysis on entire slice
print("Starting particle analysis on whole slice");
particleAnalysis(resultDir, thresholded_img, threshold_path, wholeSliceROI, closeImg, prefixwhole, filterSize);
print("Done!");

// define SCN nucleus ROI
SCN_roi_path = resultDir + File.separator + title+"_SCN_nuclei.roi";
if(File.exists(SCN_roi_path)){
	SCN_roi = SCN_roi_path;
	selectImage(originalImg);
	Stack.setDisplayMode("color");
	Stack.setChannel(channel_gcamp);
	} else {
print("Getting SCN area from user");
SCN_roi = createROI(originalImg, channel_gcamp, title, resultDir);
	}
// get Z axis profile of both whole slice and SCN nucleus
print("Extract Z-axis profile from whole slice");
extract_profile(originalImg, title, resultDir, wholeSliceROI, prefixwhole, fileName);
print("Extract Z-axis profile from SCN nuclei");
extract_profile(originalImg, title, resultDir, SCN_roi, prefixscn, fileName);

// perform aggregated particle analysis on SCN nucleus
closeImg = true;
print("Starting particle analysis on SCN nucleus");
particleAnalysis(resultDir, thresholded_img, threshold_path, SCN_roi, closeImg, prefixscn, filterSize);
print("Done!");

print("Saving Log...");
selectWindow("Log");
saveAs("Text", resultDir + File.separator + title +"_Log.txt");

// close original Image
selectImage(originalImg);
close();

print(title+ " has been analyzed!");

//anotherimage = getBoolean("Open another image");
selectWindow("Log");
run("Close");
	} // if condition end
} // for cycle end

print("Analysis completed!");
selectWindow("Log");

// FUNCTIONS


function defineSlice(channel, lowthresh, title, resultDir, originalID){ // function to identify the area of the slice
	// Args list:
		// channel		channel of interest	
		// title		title of original image
		// resultDir	path of directory where to save files
		// originalID	ID of original image
	
// identify area of slice 
	//select image
	selectImage(originalID);
	Stack.setActiveChannels(channel);
	Stack.setSlice(20);
	run("Measure");
	thresh1 = getResult("Mean", 0);
	Stack.setSlice(44);
	run("Measure");
	thresh2 = getResult("Mean", 1);
	close("Results");
	lowthresh = (thresh1+thresh2)/2;
	selectImage(originalID);
	// duplicate red channel
	run("Duplicate...", "title=masking duplicate channels=" + channel);
	
	// make binary ( title is MASK_ old title)
	//run("Make Binary", "method=Default background=Default calculate black create");
	setThreshold(lowthresh, 65535, "raw");
	setOption("BlackBackground", true);
	run("Convert to Mask", "background=Dark black");
	
	
	// expand selection
	//setOption("BlackBackground", true);
	//run("Dilate", "stack");
	
	// fill holes
	run("Fill Holes", "stack");
	
	// dilate
	run("Dilate", "stack");
	
	// analyze particles over 15000 px
	run("Analyze Particles...", "size=15000-Infinity show=Nothing clear summarize add stack");
	
	// get largest area
	selectWindow("Summary of "+ "masking"); // to edit to add real title
	
	// get all areas values
	areas = Table.getColumn("Total Area");
	max_area = Array.findMaxima(areas, 0);
	slice = max_area[0]+1;
	print("We will use Roi on slice:" +slice+ ". It has an area of "+ areas[max_area[0]]+ "px");
	
	//keep the ROI corresponding to max area
	ROItot = roiManager("count");
	roiManager("Select", max_area[0]);
	// save ROI
	extSliceROI = resultDir + File.separator + title+"_"+slice+".roi";
	roiManager("Save", extSliceROI);
		
	// Delete all ROIs
	roiManager("deselect");
	roiManager("Delete");
	
	// close MASK image
	close("masking");
	//close("MASK_masking");
	close("Summary of masking");
	
//	return ROI of interest
	return extSliceROI;
} // end of function


function particleAnalysis(resultDir, openImg, thresh_path, ROI, closeImg, prefix, filterSize) {
	// Args list:
		// resultDir		path of directory where to save files	
		// openImg			thresholded Image to perform particle analysis on
		// thresh_path		path to thresholded image saved in results dir
		// ROI				path to ROI of interest to perform analysis on
		// closeImg			true/false value. close the thresholded image?
		// prefix			prefix related to ROI being used
		
		
	// load or select image to analyse
	if (openImg != 0) {
		selectImage(openImg);
		thresh_title = getTitle();
	} else {
		open(thresh_path);
		openImg = getImageID();
		thresh_title = getTitle();
	}
	
	// load ROI of interest (select whole slice or nuclei) and clear outside area
	roiManager("open", ROI);
	roiManager("select", 0);
	selectImage(openImg);
	run("Clear Outside", "stack");
	
	// analyze particles
	run("Analyze Particles...", "size="+filterSize+"-Infinity display include clear summarize add composite stack");
	
	// save tables
	
	selectWindow("Summary of "+thresh_title);
	saveAs("Results", resultDir + File.separator + prefix + "_Summary_puncta_"+title+".csv");
	run("Close");
	
	selectWindow("Results");
	saveAs("Results", resultDir + File.separator + prefix + "_Individual_puncta_"+title+".csv");
	
	// close tables
	close("Results");
	
	// save ROIs
	roiManager("save", resultDir + File.separator + prefix + "_puncta_ROIs"+ title+"_"+".zip");
	roiManager("delete");
	if (closeImg) {
		selectImage(openImg);
		close();
		// reset global variable
		thresholded_img = -1;
		}
	
	}
	
// function for Laplace thresholding on stacks. Returns a binary image
function laplaceThresholdStack(original_img, channel, area, ROIpath, resultDir, closeImg, threshold) {
	
	// Args list:
		// original_img			ID of unprocessed 3 channel img
		// channel				channel to be duplicated
		// area					true/false is there an area to be analysed
		// ROIpath				path to ROI to be loaded
		// resultDir			path to directory where to store results
		// closeImg				true/false do you want to close the thresholded img at the end
		// threshold			value to use for thresholding image
		
	// select imageID and duplicate channel of interest
    selectImage(original_img);
	original_title = File.nameWithoutExtension();
	run("Duplicate...", "title=to_process duplicate channels="+ channel);
	to_process = getImageID();
	
	//check if user has selected the an area to restrict the analysis to.
	if (area){
		roiManager("open", ROIpath);
		}else{
			run("Select All");
			roiManager("add");
			}
	
	// clear outside
	roiManager("select", 0);
	run("Clear Outside", "stack");
	
    selectImage(to_process);
    // Get the number of slices in the stack
    stackSize = nSlices+1;
    setBatchMode(true);
    // Create a new stack for the Laplace-transformed images
    newImage("Laplace_Stack", "8-bit black", getWidth(), getHeight(), 1);
    selectWindow("Laplace_Stack");
    laplaceStackTitle = getTitle();
    
    // Loop through each slice in the stack
    for (i = 1; i < stackSize; i++) {
        selectImage(to_process); // Select the original stack
        setSlice(i);                 // Go to the current slice
        
        // Duplicate the current slice for processing
        run("Duplicate...", "title=Working_Slice");
        selectWindow("Working_Slice");
        
        // Apply Laplacian kernel
        run("Convolve...", "text1=[-1 -1 -1 -1 -1\n-1 1 1 1 -1\n-1 1 8 1 -1\n-1 1 1 1 -1\n-1 -1 -1 -1 -1\n] normalize"); // Laplacian kernel
        
        // Get statistics for thresholding
        getStatistics(area, mean, min, max, stdDev);
        thresholdValue = -threshold * stdDev; // Normalize threshold by negative SD
        
        // Apply threshold
        setAutoThreshold("Default dark");
        setThreshold(threshold, max, "dark");
        run("Convert to Mask");
        
        // Save Laplacian-transformed slice to the Laplace stack
        run("Concatenate...", "  title=Laplace_Stack open image1=Laplace_Stack image2=Working_Slice image3=[-- None --]");
        
        // Cleanup intermediate windows
        close("Working_Slice");
    }
    selectWindow(laplaceStackTitle);
    setSlice(1);
    run("Delete Slice");
    setBatchMode(false);
    
    // save image
    selectWindow(laplaceStackTitle);
	thresholded_img = getImageID();
    threshold_path = resultDir + File.separator + original_title + "_thresholded.tiff";
	saveAs("Tiff", threshold_path);
	
	if (closeImg) {
		close();
		// reset global variable
		thresholded_img = 0;
		}
    // Return to the original stack
    selectImage(to_process);
    close();
    return threshold_path;
}

function threshold_prep(original_img, channel, area, ROIpath, resultDir, gammasave, closeImg){ // function to prepare the image for particle analysis
	// Args list:
		// original_img			ID of unprocessed 3 channel img
		// channel				channel to be duplicated
		// area					true/false is there an area to be analysed
		// ROIpath				path to ROI to be loaded
		// resultDir			path to directory where to store results
		// gammasave			true/false do you want to save the image before thresholding?
		// closeImg				true/false do you want to close the thresholded img at the end
		
	// select imageID and duplicate channel of interest
	selectImage(original_img);
	original_title = File.nameWithoutExtension();
	run("Duplicate...", "title=to_process duplicate channels="+ channel);
	to_process = getImageID();
	
	//check if user has selected the an area to restrict the analysis to.
	// this allows to get the median values from the wholw image or the slice
	// area only. The particle analysis will be done subsequently on a
	// different ROI selected by the code
	if (area){
		roiManager("open", ROIpath);
		}else{
			run("Select All");
			roiManager("add");
			}
	
	// clear outside
	roiManager("select", 0);
	run("Clear Outside", "stack");
	
	// get slices count
	nslices = nSlices;
	
	// get measurements from all slices
	for (i = 0; i < nslices; i++) {
	roiManager("select", 0);
	j = i+1;
	setSlice(j);
	run("Measure");
	// save first row of table to access all vals
	if (i == 0) {
		selectWindow("Results");
		saveAs("Results", resultDir + File.separator + original_title + "_slice_params.csv");
		}
	}
	
	
	
	// extract median and StdDev values
	median = Table.getColumn("Median");
	stdDev = Table.getColumn("StdDev");
	mean = Table.getColumn("Mean");
	max = Table.getColumn("Max");
	
	// subtract values of combined median and StdDev values from image
	sub_array = Array.getSequence(nslices);
	for (i = 0; i < nslices; i++) {
	j = i+1;
	setSlice(j);
	subtracted = median[i] - stdDev[i];
	run("Subtract...", "value="+ subtracted +" slice");
	
	// add values to array
	sub_array[i] = subtracted;
	}
	
	// create new table with processed values
	Table.create("Sub_values");
	// add columns
	Table.setColumn("Mean", mean);
	Table.setColumn("Max", max);
	Table.setColumn("Median", median);
	Table.setColumn("StdDev", stdDev);
	Table.setColumn("Subtracted", sub_array);
	
	// save table and close
	saveAs("Results", resultDir + File.separator + original_title + "_sub_values.csv");
	close(original_title + "_sub_values.csv");
	
	// Apply gamma and Threshold image
	selectImage(to_process);
	run("Gamma...", "value=1.50 stack");
	if (gammasave) {saveAs("Tiff", resultDir + File.separator + original_title + "_processed.tiff");}
	setThreshold(60, 65535, "raw");
	thresholded_img = getImageID();
	
	// save image and close (optional)
	threshold_path = resultDir + File.separator + original_title + "_thresholded.tiff";
	saveAs("Tiff", threshold_path);
	
	if (closeImg) {
		close();
		// reset global variable
		thresholded_img = 0;
		}
	// close ROI manager
	roiManager("delete");
	selectWindow("ROI Manager");
	run("Close");
	
	return threshold_path
}


function createROI(originalImg, channel, title, resultDir){ // function to create a good ROI around nucleus
	// Args list:
		// originalImg	ID of the original image
		// channel		channel to display
		// title		title of original image
		// resultDir	path to directory where to store results
	
	// select original image and activate channel
	selectImage(originalImg);
	Stack.setDisplayMode("color");
	Stack.setChannel(channel);
	
	// while cycle to get a nice ROI
	niceROI = false;
	while (!niceROI) {
	if (roiManager("count") != 0) {
		roiManager("delete");
	}
	
	// draw ROI and press OK
	setTool("freehand");
	waitForUser("Draw ROI around Nuclei and press OK.");
	roiManager("add");
	
	// ask if the ROI is good enough
	niceROI = getBoolean("Are you happy with the ROI selection?");
	}
	
	// save ROI
	roiManager("Select", 0);
	saveROIpath = resultDir + File.separator + title+"_SCN_nuclei.roi";
	roiManager("Save", saveROIpath);
	
	// clean ROImanager and return path to ROI
	roiManager("delete");
	return saveROIpath;
} //end of createROI function

// write foo to extract and save Zaxis profiles
function extract_profile(originalImg, title, resultDir, ROI, prefix, filename){
	// Args list:
		// originalImg	ID of original image
		// title		title of original image
		// resultDir	path to directory where to store results
		// ROI			path to ROI to use for Z-axis profile
		// prefix		prefix to use in file saving
	
	// Select image
	selectImage(originalImg);
	// activate all channels
	Stack.setDisplayMode("composite");
	// open ROI of interest
	roiManager("open", ROI);
	// get multichannel ZT axis profile
	run("Multichannel ZT-axis Profile", "statschoice=Mean calibratedx=true zoption=[Active plane] activatechannelwidget=true allowroutine=true");
	plotWin = "Z/T-axis Plot of "+filename;
	selectWindow(plotWin);
	// create table with values
	Plot.showValuesWithLabels();
	// change column labels
	Table.renameColumn("Mean_(Active_plane)", "FRED");
	Table.renameColumn("Y1", "GCaMP8s");
	//Table.renameColumn("Y2", "phase");
	// save values and close table
	saveAs("Results", resultDir + File.separator + prefix + "_Z_profile_multi_"+title+".csv");
	close("Results");
	close(plotWin);
	
	// close ROI manager
	selectWindow("ROI Manager");
	run("Close");
	}

function setImgProperties(pxSize){
//	args
	// pxSize	pixelSize 
	
	Stack.setXUnit("um");
	run("Properties...", "pixel_width="+pxSize+" pixel_height="+pxSize+" voxel_depth="+pxSize+" frame=[1800 sec]");
	}

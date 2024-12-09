// macro to extract crop of SCN slices (Left/right) 

//get folder containing files
inputFold = getDir("Select directory where files are stored");
// get list of files
filelist = getFileList(inputFold) 
// ask if you need to create the Directoriess
createDirs = "N";
resaveROIs = "Y";
adjustROI = "N";
createDirs = getString("Do you want to create directories", "Y");

// ask if you need to create the ROIs
createROIs = "Y";
createROIs = getString("Do you need to define ROIs", "Y");
if (createROIs != "Y") {
	adjustROI = getString("Do you need to adjust ROIs", "Y");
	}
// ask if you need to extract SCNs
extractSCNs = "Y";
extractSCNs = getString("Do you need to extract SCN stacks", "Y");
extractNuclei = getString("Do you want to extract SCN nuclei only?", "N");

// start loop for processing files
for (i = 0; i < lengthOf(filelist); i++) {
		    if (endsWith(filelist[i], ".tif") | endsWith(filelist[i], ".tiff")) {
		    	filepath = inputFold + File.separator + filelist[i];
		        open(filepath);
		        
		// get filename without extension
		filename = File.nameWithoutExtension();
		fullname = File.getName(filepath);
		// get fileID
		fileID = getImageID();
		
		if (createDirs == "Y") {
		
			// create folder of results 
			newDir = inputFold + File.separator + "hemis";
			dirExist = File.isDirectory(newDir);
			if (dirExist != 1) {
				print("Creating directory to store images");
				File.makeDirectory(newDir);
			}
			}
		
		//if create ROIs = "Y"
// this section is to create ROIs for the SCN
		if (createROIs == "Y") {
		
			// create folder of results 
			newDir = inputFold + File.separator + "hemis";
			dirExist = File.isDirectory(newDir);
			if (dirExist != 1) {
				print("Creating result directory");
				File.makeDirectory(newDir);
			}
			// empty ROI manager

			roiManager("reset");
			Stack.setActiveChannels("010");
			Stack.setChannel(2)
			run("Enhance Contrast", "saturated=0.35");
			// ask user to draw ROI around left side of slice
			setTool("rectangle");
		    waitForUser("Draw an ROI around the left side of the slice and press ok.");
			// add ROI and rename as LEFT_SLICE_filename
			roiManager("add");
			name = "LEFT_"+filename;
			roiManager("select", 0);
			roiManager("rename", name);
			roiManager("add");
			roiManager("select", 1);
			// ask user to draw ROI around right side of slice

			setTool("rectangle");
		    waitForUser("Draw an ROI around the right side of the slice and press ok.");
		    // add ROI and rename as RIGHT_SLICE_filename
			roiManager("update");
			name2 = "RIGHT_"+filename;
			roiManager("select", 1);
			roiManager("rename", name2);
			// save ROI as zip within result folder
			ROIzippath =  inputFold + File.separator + "hemis" + File.separator + filename+"_ROIset.zip";
			roiManager("save", ROIzippath);
		}else {
			// open the archive containing the ROIs
			ROIpath = inputFold + File.separator + "hemis" + File.separator + filename + "_ROIset.zip";
			roiManager("reset");
			// check if the ROI file exists at the created path otherwise ask the user to select the ROI file
			if(File.exists(ROIpath)){
			roiManager("open", ROIpath);
			}else {
			Stack.setActiveChannels("010");
			Stack.setChannel(2)
			zipFile = File.openDialog("Choose a ROIset.zip file");
			open(zipFile);
			if (adjustROI == "Y") {
				roiManager("select", 0);
				waitForUser("Adjust left ROI and press ok.");
				roiManager("update");
				roiManager("select", 1);
				waitForUser("Adjust right ROI and press ok.");
				roiManager("update");
				}
			if (resaveROIs == "Y") {
			// save ROI as zip within result folder
			ROIzippath =  inputFold + File.separator +"hemis" + File.separator + filename+"_ROIset.zip";
			roiManager("save", ROIzippath);
			}
			}
		}
		// if extractSCNs = "Y"
// this section is to create the files isolating only specific areas of the slice
		if (extractSCNs == "Y") {
			LeftDir = newDir + File.separator + "Left";
			dirExist = File.isDirectory(LeftDir);
			if (dirExist != 1) {
				print("Creating result directory");
				File.makeDirectory(LeftDir);
			}
			
			// select image
			selectImage(fileID);
			
			// create left slice image: select ROI, duplicate stack, clear outside, flip, select new ROI
			selectImage(fileID);
			roiManager("select", 0);
			run("Duplicate...", "title=LEFT_"+fullname+" duplicate channels=1-2");
			setBackgroundColor(0, 0, 0);
			run("Flip Horizontally");
			if (extractNuclei == "Y") {
				if (adjustROI == "Y") {
					roiManager("select", 2);
				waitForUser("Adjust ROI and press ok.");
				roiManager("update");
					}else {
			setTool("freehand");
		    waitForUser("Draw an ROI around the nucleus and press ok.");
		    roiManager("add");
					}
		    run("Clear Outside", "stack");
		    roiManager("select", 2);
		    name3 = "LEFT_SCN_"+filename;
		    roiManager("rename", name3);
			}
			SliceStackPath = LeftDir + File.separator + "LEFT_"+fullname;
			saveAs("Tiff", SliceStackPath);
			close("*LEFT*");
						
		// create right slice image: select ROI, duplicate stack, clear outside, save
		
		// create right directory
		RightDir = newDir + File.separator + "Right";
			dirExist = File.isDirectory(RightDir);
			if (dirExist != 1) {
				print("Creating result directory");
				File.makeDirectory(RightDir);
			}
		
			selectImage(fileID);
			roiManager("select", 1);
			run("Duplicate...", "title=RIGHT_"+fullname+" duplicate channels=1-2");
			if (extractNuclei == "Y") {
				if (adjustROI == "Y") {
					roiManager("select", 3);
				waitForUser("Adjust ROI and press ok.");
				roiManager("update");
					}else {
			setTool("freehand");
		    waitForUser("Draw an ROI around the nucleus and press ok.");
		    roiManager("add");
					}
		    run("Clear Outside", "stack");
		    roiManager("select", 3);
		    name4 = "RIGHT_SCN_"+filename;
		    roiManager("rename", name4);
		    }
			setBackgroundColor(0, 0, 0);
			
			SliceStackPath = RightDir + File.separator + "RIGHT_"+fullname;
			saveAs("Tiff", SliceStackPath);
			close("*RIGHT*");
			
			//Re-Save ROIs
			roiManager("Deselect");
			ROIzippath =  inputFold + File.separator +"hemis" + File.separator + filename+"_ROIset.zip";
			roiManager("save", ROIzippath);
		} // end of if extractSCN =="Y"
		selectImage(fileID);
		close();
		
		} //end of if is tiff

} // end of loop for opening files
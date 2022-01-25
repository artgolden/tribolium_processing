//Stack Projection
setBatchMode(true);

inputDir = getDirectory( "Choose the Input Directory" );

dataset_name = getString("Enter the dataset identifier (will be prepended to the output file name).");

top_embryo_direction1_ID = "SPC00" + getString("ID corresponding to the TOP embryo direction 1 (SPC00N), enter N:", "0");
top_embryo_direction2_ID = "SPC00" + getString("ID corresponding to the TOP embryo direction 2 (SPC00N), enter N:", "1");
top_embryo_direction3_ID = "SPC00" + getString("ID corresponding to the TOP embryo direction 3 (SPC00N), enter N:", "2");
top_embryo_direction4_ID = "SPC00" + getString("ID corresponding to the TOP embryo direction 4 (SPC00N), enter N:", "3");
down_embryo_direction1_ID = "SPC00" + getString("ID corresponding to the DOWN embryo direction 1 (SPC00N), enter N:", "4");
down_embryo_direction2_ID = "SPC00" + getString("ID corresponding to the DOWN embryo direction 2 (SPC00N), enter N:", "5");
down_embryo_direction3_ID = "SPC00" + getString("ID corresponding to the DOWN embryo direction 3 (SPC00N), enter N:", "6");
down_embryo_direction4_ID = "SPC00" + getString("ID corresponding to the DOWN embryo direction 4 (SPC00N), enter N:", "7");

list=getFileList(inputDir); 
Array.sort(list);

for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], "TIF") || endsWith(list[i], "tif")) {
		print(list[i]);
		open(inputDir+list[i]);
		//run("Bio-Formats Importer", "open=[inputDir+list[i]] color_mode=Default view=[Standard ImageJ] stack_order=Default");
		fName=File.nameWithoutExtension;
		print(fName);
		input=getTitle();
		Stack.getDimensions(imageWidth, imageHeight, numChannels, numSlices, frames);
		resize = 4 * numSlices;
		
		
		// Top embryo direction 1 
		if(matches(list[i], "(.*)" + top_embryo_direction1_ID + "(.*)(TIF|tif)")){
			direction = "DR0001"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo");
			File.makeDirectory(outputDir);
			
		}
		// Top embryo direction 2
		if(matches(list[i], "(.*)" + top_embryo_direction2_ID + "(.*)(TIF|tif)")){ 
			direction = "DR0002"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo");
			File.makeDirectory(outputDir);
		}
		// Top embryo direction 3
		if(matches(list[i], "(.*)" + top_embryo_direction3_ID + "(.*)(TIF|tif)")){ 
			direction = "DR0003"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo");
			File.makeDirectory(outputDir);
		}
		// Top embryo direction 4
		if(matches(list[i], "(.*)" + top_embryo_direction4_ID + "(.*)(TIF|tif)")){ 
			direction = "DR0004"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "top_embryo");
			File.makeDirectory(outputDir);
		}
		// Down embryo direction 1 
		if(matches(list[i], "(.*)" + down_embryo_direction1_ID + "(.*)(TIF|tif)")){
			direction = "DR0001"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo");
			File.makeDirectory(outputDir);
		}
		// Down embryo direction 2
		if(matches(list[i], "(.*)" + down_embryo_direction2_ID + "(.*)(TIF|tif)")){ 
			direction = "DR0002"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo");
			File.makeDirectory(outputDir);
		}
		// Down embryo direction 3
		if(matches(list[i], "(.*)" + down_embryo_direction3_ID + "(.*)(TIF|tif)")){ 
			direction = "DR0003"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo");
			File.makeDirectory(outputDir);
		}
		// Down embryo direction 4
		if(matches(list[i], "(.*)" + down_embryo_direction4_ID + "(.*)(TIF|tif)")){ 
			direction = "DR0004"
			outputDir = File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo" + File.separator + direction + File.separator;
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp");
			File.makeDirectory(File.getParent(inputDir) + File.separator + "mp" + File.separator + "down_embryo");
			File.makeDirectory(outputDir);
		}

		file_index = parseInt(substring(fName, indexOf(fName, "TL") + 2, indexOf(fName, "TL") + 6)) + 1;
		new_file_name = dataset_name + "TP" + IJ.pad(file_index, 4) + direction + "CH0001PL(ZM)"
	

		
		//run("Properties...", "channels=numChannels slices=numSlices frames=1 unit=µm pixel_width=0.645 pixel_height=0.645 voxel_depth=2.58 frame=[0 sec] origin=0,0");
		run("Z Project...", "start=1 stop=numSlices projection=[Max Intensity]");
		saveAs("Tiff", outputDir + new_file_name + "XYMAX.tif");
		//run("Scale Bar...", "width=100 height=12 font=44 color=White background=None location=[Lower Right]");
		//saveAs("Tiff", outputDir+fName+"XYMAXScale.tif");
	
		//selectWindow(input);
		//run("Reslice [/]...", "output=2.580 start=Bottom avoid");
		//run("Z Project...", "start=1 stop=imageHeight projection=[Max Intensity]");
		//run("Scale...", "x=1.0 y=4 width="+imageWidth+" height="+imageHeight+" interpolation=None average create title=[xz.tif]");
		//saveAs("Tiff", outputDir+fName+"XZMAX.tif");
		//run("Scale Bar...", "width=100 height=12 font=44 color=White background=None location=[Lower Right]");
		//saveAs("Tiff", outputDir+fName+"XZMAXScale.tif");

		//selectWindow(input);
		//run("Reslice [/]...", "output=2.580 start=Right avoid");
		//run("Z Project...", "start=1 stop=imageWidth projection=[Max Intensity]");
		//run("Scale...", "x=1.0 y=4 width="+imageWidth+" height="+imageHeight+" interpolation=None average create title=[zy.tif]");
		//saveAs("Tiff", outputDir+fName+"ZYMAX.tif");
		//run("Scale Bar...", "width=100 height=12 font=44 color=White background=None location=[Lower Right]");
		//saveAs("Tiff", outputDir+fName+"ZYMAXScale.tif");

		while (nImages>0) {
			selectImage(nImages); 
			close(); 
		}	
	}
} 

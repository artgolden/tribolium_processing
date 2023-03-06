#@ File (label="Select a directory", style="directory") myDir
#@ File (label="Select a output directory", style="directory") outputDir
#@ Integer (label="Number of channels", min=0, max=10, value=3) num_channels
#@ String (label="Channels substitution string", description="Name field") ch_to_keep

from __future__ import print_function
from ij import IJ, WindowManager


def main():
	root = myDir.getPath() # get the root out the java file object
	
	import os, glob

	# set up bioformats
	from loci.plugins import BF
	from loci.plugins.in import ImporterOptions
	options = ImporterOptions()
	options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
	options.setGroupFiles(True)  # if the files are named logically if will group them into a stack
  	

	for path, subdirs, files in os.walk(root):
		# just get the one of the files that matches your image pattern
		flist = glob.glob(os.path.join(path,"*CH*.tif"))
		print(flist, "  -  ", path)
		if( flist ):
			file = flist[0]
			print("Processing {}".format(file))
			options.setId(file)
			imp = BF.openImagePlus(options)[0]
			print(imp.getClass())
			if num_channels > 1:
				imp = IJ.run(imp, "Arrange Channels...", "new=%s" % ch_to_keep);
			


			IJ.save(imp, os.path.join(outputDir.getPath(), os.path.basename(file).rsplit('_',2)[0] + "_hyperstack.tif"))

			# closes the windows if they are open
			imp = WindowManager.getCurrentImage()
			imp.close()
main()
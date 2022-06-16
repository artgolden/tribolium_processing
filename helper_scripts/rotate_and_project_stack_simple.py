#@ Integer middle_plane
# @ File(label='Choose a file with a z-stack to rotate and project', style='file') imp_file
# @ File(label='Choose a directory to save projections to', style='directory') projections_dir

import os
from ij import IJ
from ij.plugin import ZProjector
from ij.plugin import Slicer
from ij.io import FileSaver

def save_tiff(image, path):
	if os.path.exists(path):
		os.remove(path)
	fs = FileSaver(image)
	fs.saveAsTiff(path)


imp_path = imp_file.getAbsolutePath()
projection_path_no_extension = os.path.join(projections_dir.getAbsolutePath(), os.path.splitext(os.path.basename(imp_path))[0])
print(projection_path_no_extension)
imp = IJ.openImage(imp_path)
projection = ZProjector.run(imp,"max",1,middle_plane);
save_tiff(projection, projection_path_no_extension + "_angle_0.tiff")
for i in range(15):
	resliced = Slicer.reslice(Slicer(), imp)
	imp.close()
	IJ.run(resliced, "Rotate... ", "angle=22.5 grid=1 interpolation=Bilinear stack");
	imp = Slicer.reslice(Slicer(), resliced)
	projection = ZProjector.run(imp,"max",1,middle_plane);
	save_tiff(projection, projection_path_no_extension + "_angle_%s.tiff" % (round(22.5 * (i + 1))))
	


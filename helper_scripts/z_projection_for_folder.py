#@ File (label='Choose image') image_path
#@ Integer (label='Start', value=1) start
#@ Integer (label='End', value=489) end
#@ File (label='Choose an output directory', style='directory') output_dir



''' A script to make Z-projection for a stack. Applicable to a folder with batch processing.'''

import os
from ij import IJ
from ij.io import FileSaver
from ij.plugin import ZProjector

def save_tiff_simple(image, path):
    if os.path.exists(path):
        os.remove(path)
    fs = FileSaver(image)
    fs.saveAsTiff(path)
image_path = image_path.getAbsolutePath()
output_file_path = os.path.join(output_dir.getAbsolutePath(), os.path.basename(image_path))


imp = IJ.openImage(image_path)
imp = ZProjector.run(imp,"max",start,end)
save_tiff_simple(imp, output_file_path)
#@ File (label='Choose image to expand canvas') image_path
#@ File (label='Choose an output directory', style='directory') output_dir


''' A script to make first plane in a z-stack black. Applicable to a folder with batch processing.'''

import os
from ij import IJ
from ij.io import FileSaver

def save_tiff_simple(image, path):
    if os.path.exists(path):
        os.remove(path)
    fs = FileSaver(image)
    fs.saveAsTiff(path)
image_path = image_path.getAbsolutePath()
output_file_path = os.path.join(output_dir.getAbsolutePath(), os.path.basename(image_path))


imp = IJ.openImage(image_path)
IJ.run(imp, "Select All", "")
IJ.setForegroundColor(0, 0, 0)
IJ.run(imp, "Fill", "slice")
save_tiff_simple(imp, output_file_path)
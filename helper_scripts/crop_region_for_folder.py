#@ File (label='Choose image to expand canvas') image_path
#@ File (label='Choose an output directory', style='directory') output_dir


''' A script to crop region from a z-stack. Applicable to a folder with batch processing.'''

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
imp.setRoi(94,311,450,240)
imp = imp.resize(450, 240, 1, "bilinear")
imp = imp.crop("35-95")
save_tiff_simple(imp, output_file_path)
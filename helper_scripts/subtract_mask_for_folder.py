#@ File (label='Choose image') image_path
#@ File (label='Choose mask to use for subtraction') mask_path
#@ File (label='Choose an output directory', style='directory') output_dir
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert

''' A script to make subtract a mask from a z-stack. Applicable to a folder with batch processing.'''

import os
from ij import IJ
from ij.io import FileSaver
from ij.plugin import ImageCalculator

def save_tiff_simple(image, path):
    if os.path.exists(path):
        os.remove(path)
    fs = FileSaver(image)
    fs.saveAsTiff(path)
image_path = image_path.getAbsolutePath()
mask_path = mask_path.getAbsolutePath()
output_file_path = os.path.join(output_dir.getAbsolutePath(), os.path.basename(image_path))


imp = IJ.openImage(image_path)
mask = IJ.openImage(mask_path)
IJ.run(mask, "16-bit", "");
IJ.run(mask, "Multiply...", "value=256 stack")
subtracted = ImageCalculator.run(imp, mask, "Subtract stack")
save_tiff_simple(subtracted, output_file_path)
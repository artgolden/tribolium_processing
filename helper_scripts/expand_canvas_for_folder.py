#@ File (label='Choose image to expand canvas') image_path
#@ File (label='Choose an output directory', style='directory') output_dir
#@ Integer (label='How much to expand width (px)') width_expansion
#@String(label='How to center image', choices={"Center-Left", "Center-Right"}, style="listBox", value="Center-Right") image_centering_pos

''' A script to expand canvas width of the image. Applicable to a folder with batch processing.'''

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
IJ.run(imp, "Canvas Size...", "width=%s height=%s position=%s zero" % (imp.getWidth() + width_expansion, imp.getHeight(), image_centering_pos))
save_tiff_simple(imp, output_file_path)
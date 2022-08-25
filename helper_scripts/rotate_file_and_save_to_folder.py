#@ File image_to_rotate
#@ File(label='Directory to save images', style='directory') output_dir
#@ String (label='How to rotate the image?', choices={"Right", "Left", "Flip Vertically", "Flip Horizontally"}, style="radioButtonHorizontal", value="Left") rotation_direction

image_to_rotate = image_to_rotate
output_dir = output_dir
rotation_direction = rotation_direction

import os
from ij import IJ
from ij.io import FileSaver

def save_tiff_simple(image, path):
    if os.path.exists(path):
        os.remove(path)
    fs = FileSaver(image)
    fs.saveAsTiff(path)

image_path = image_to_rotate.getAbsolutePath()
output_path = os.path.join(output_dir.getAbsolutePath(), os.path.basename(image_path))
print("Opening image %s" % image_path)
imp = IJ.openImage(image_path)
if rotation_direction in ["Flip Vertically", "Flip Horizontally"]:
    IJ.run(imp, rotation_direction, "stack")
else:
    IJ.run(imp, "Rotate 90 Degrees %s" % rotation_direction, "")
print("Saving rotated to the %s image %s" % (rotation_direction, output_path))
save_tiff_simple(imp, output_path)



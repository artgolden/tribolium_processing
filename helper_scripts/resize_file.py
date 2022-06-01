# @ File(label='File to resize', style='file') imp_path
# @ Integer (label='dummy variable', value=1) dummy


from ij import IJ
from ij.io import FileSaver

imp_path = imp_path.getAbsolutePath()
imp = IJ.openImage(imp_path)
imp = imp.resize(1800, 1200, "bilinear")

fs = FileSaver(imp)
fs.saveAsTiff(imp_path)
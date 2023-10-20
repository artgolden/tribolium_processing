#@ ImagePlus imp
#@ Integer (label='Scale factor', value=8) scale_factor


''' A script to make a reslice from top and then Z-projection for a stack. Used for bead images quality checking.'''

import os
from ij import IJ, WindowManager
from ij.plugin import ZProjector


num_planes = imp.getNSlices()
reslice = IJ.run(imp, "Reslice [/]...", "output=1.000 start=Top");
reslice = WindowManager.getCurrentImage()
reslice.hide()
reslice = ZProjector.run(reslice,"max");
print num_planes
reslice = reslice.resize(imp.getWidth(), num_planes * scale_factor, "bilinear");
reslice.show()
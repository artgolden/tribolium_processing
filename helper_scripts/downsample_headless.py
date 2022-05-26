#@ ImagePlus imp
#@ Integer (label='width') width
#@ Float (label='sourceSigma', value=0.5) sourceSigma
#@ Float (label='targetSigma', value=0.5) targetSigma
#@output ImagePlus downsampled_imp

import math
from ij import IJ


def downsample_image_stack(stack, new_width, sourceSigma=0.5, targetSigma=0.5):
    """Downsample image stack given its width

    Args:
        stack (ImagePlus): image stack to downsample
        new_width (int): new width
        sourceSigma (float, optional): source image sigma for gaussian blur. Defaults to 0.5.
        targetSigma (float, optional): target image sigma for gaussian blur. Defaults to 0.5.

    Returns:
        ImagePlus: downsampled image stack
    """
    stack = stack.duplicate()
    height = int(round(float(new_width) * stack.getHeight() / stack.getWidth()))

    if new_width <= stack.getWidth():
        s = targetSigma * stack.getWidth() / new_width

        IJ.run(stack, "Gaussian Blur...", "sigma=%s stack" % math.sqrt(s * s - sourceSigma * sourceSigma))
        IJ.run(stack, "Scale...", "x=- y=- width=%s height=%s process title=- interpolation=None" % (new_width, height))

        extraX = 0 if stack.getWidth() % 2 == 0 else 1
        extraY = 0 if stack.getHeight() % 2 == 0 else 1
        initialX = (stack.getWidth() / 2 - new_width / 2 + extraX) if new_width % 2 == 0 else (stack.getWidth() / 2 - new_width / 2 + 1 - extraX)
        initialY = (stack.getHeight() / 2 - height / 2 + extraY) if height % 2 == 0 else (stack.getHeight() / 2 - height / 2 + 1 - extraY)
        stack.setRoi(initialX, initialY, new_width, height) 
        downsampled_imp = stack.crop("stack")
    else:
        IJ.log("You try to upsample the image.  You need an interpolator for that not a downsampler.")
    return downsampled_imp


downsampled_imp = downsample_image_stack(imp, width, sourceSigma, targetSigma)

#@ ImagePlus mask



from distutils.dir_util import mkpath
import math
import os
import re
import fnmatch
import json
import logging
import traceback


from java.io import File
from ij.io import FileSaver
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin.filter import RankFilters
from fiji.threshold import Auto_Threshold
from ij.plugin.filter import ParticleAnalyzer
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from ij.plugin import RoiRotator
from ij.plugin import ZProjector
from ij.plugin import Slicer
from ij.plugin import StackCombiner
from ij.plugin import StackMaker
from ij.process import ByteProcessor
from ij.io import RoiDecoder
from ij.gui import PointRoi, RotatedRectRoi
from emblcmci import BleachCorrection_MH

def smoothen_embryo_mask(mask, num_smoothing_iterations):
    # Repeated Erosion and Dilation to remove dirt, which shows as protrusions on the embryo mask
    # Using IJ.run() for erosion and dilation has some persistent state after run
    # that produces inconsistent behaviour and bugs
    # This is why I am using ByteProcessor methods directly
    bt = mask.getProcessor().convertToByteProcessor()
    
    for i in range(num_smoothing_iterations):
        bt.erode(2, 0)
    for i in range(num_smoothing_iterations):
        bt.dilate(2, 0)
    for i in range(num_smoothing_iterations):
        bt.erode(3, 0)
    for i in range(num_smoothing_iterations):
        bt.dilate(3, 0)
    mask.setProcessor(bt)


def bounding_roi_and_params_from_embryo_mask(mask):
    IJ.run("Clear Results")
    table = ResultsTable()
    # Initialise without display (boolean value is ignored)
    roim = RoiManager(False)
    MIN_PARTICLE_SIZE = (int(mask.getWidth()) / 3) ** 2  # pixel ^ 2
    MAX_PARTICLE_SIZE = float("inf")
    ParticleAnalyzer.setRoiManager(roim)
    pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER,
            Measurements.ELLIPSE,
            table,
            MIN_PARTICLE_SIZE,
            MAX_PARTICLE_SIZE)
    pa.analyze(mask)

    rot_angle = round(table.getValue("Angle", 0), ndigits=1)
    if rot_angle > 90:
        rot_angle -= 180
    print("\t Angle of the elipse fitted onto the embryo: %s", rot_angle)

    roi_arr = roim.getRoisAsArray()
    roim.runCommand('reset')
    mask.setRoi(roi_arr[len(roi_arr) - 1])
    IJ.run(mask, "Select Bounding Box (guess background color)", "")
    bounding_roi = mask.getRoi()
    # Crop the image so the borders are equidistant from the borders of the bounding box
    # for the center of rotation to be aligned roughly with the center of the embryo
    bounding_rectangle = bounding_roi.getBounds()
    # extracting width field from java.awt.Rectangle
    i = bounding_rectangle.x
    j = bounding_rectangle.y
    We = bounding_rectangle.width
    He = bounding_rectangle.height
    W = mask.getWidth()
    H = mask.getHeight()
    ax = i
    bx = W - i - We
    ay = j
    by = H - j - He
    xn = i - min(ax, bx)
    yn = j - min(ay, by)
    Wn = We + 2 * min(ax, bx)
    Hn = He + 2 * min(ay, by)
    IJ.run(mask, "Specify...", "width=%s height=%s x=%s y=%s" % (Wn, Hn, xn, yn))
    equidistant_crop_roi = mask.getRoi()
    mask.setRoi(equidistant_crop_roi)
    mask_cropped = mask.crop()
    # For some reason cropping a mask does not produce a proper binarized image (it is gray if you view it), 
    # So I had to redo the thresholding
    IJ.run(mask_cropped, "Rotate... ", "angle=%s grid=1 interpolation=Bilinear" % rot_angle)
    hist = mask_cropped.getProcessor().getHistogram()
    triag_threshold = Auto_Threshold.Triangle(hist)
    mask_cropped.getProcessor().setThreshold(
        triag_threshold, float("inf"), ImageProcessor.NO_LUT_UPDATE)
    IJ.run(mask_cropped, "Convert to Mask", "")
    roim.close()

    # Extract the precise embryo center in the rotated cropped image
    roim = RoiManager(False)
    MIN_PARTICLE_SIZE = (int(mask.getWidth()) / 3) ** 2  # pixel ^ 2
    MAX_PARTICLE_SIZE = float("inf")
    ParticleAnalyzer.setRoiManager(roim)
    # For some reason ParticleAnalyzer cannot be run without reinitialization,
    # that's why repeating it here
    pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER,
            Measurements.ELLIPSE,
            table,
            MIN_PARTICLE_SIZE,
            MAX_PARTICLE_SIZE)
    pa.analyze(mask_cropped)
    table.reset()
    roi_arr = roim.getRoisAsArray()
    mask_cropped.setRoi(roi_arr[len(roi_arr) - 1])
    roim.runCommand('reset')
    roim.close()

    IJ.run(mask_cropped, "Select Bounding Box (guess background color)", "")
    cropped_bounding_roi = mask_cropped.getRoi()
    embryo_length = cropped_bounding_roi.getBounds().width
    embryo_width = cropped_bounding_roi.getBounds().height
    print("\tEmbryo dims: (%s, %s)" % (embryo_length, embryo_width))
    embryo_cropped_rot_center_x = cropped_bounding_roi.getContourCentroid()[0]
    embryo_cropped_rot_center_y = cropped_bounding_roi.getContourCentroid()[1]
    # Transform the embryo center coordinates with the inverse rotation
    rad_angle = math.radians(-rot_angle)
    # switch to coordinates with origin in the cropped image center.
    xC = embryo_cropped_rot_center_x - Wn / 2
    yC = embryo_cropped_rot_center_y - Hn / 2
    embryo_cropped_center_x = xC * \
        math.cos(rad_angle) - yC * math.sin(rad_angle)
    embryo_cropped_center_y = xC * \
        math.sin(rad_angle) + yC * math.cos(rad_angle)
    embryo_cropped_center_x += Wn / 2
    embryo_cropped_center_y += Hn / 2
    # Transform the embryo center coordinates to reverse cropping
    embryo_center_x = embryo_cropped_center_x + xn
    embryo_center_y = embryo_cropped_center_y + yn
    
    
    print("\tEmbryo center: (%s, %s)" % (embryo_center_x, embryo_center_y))
    box_width = embryo_length + 12 + 4 - embryo_length % 4
    box_height = embryo_width + 12 + 4 - embryo_width % 4
    IJ.run(mask, "Specify...", "width=%s height=%s x=%s y=%s centered" %
            (box_width, box_height, embryo_center_x, embryo_center_y))
    bounding_roi = mask.getRoi()
    bounding_roi_rot = RoiRotator.rotate(bounding_roi, -rot_angle)
    bounding_roi = {
        "bounding_roi_rect" : bounding_roi_rot,
        "embryo_length" : embryo_length,
        "embryo_width" : embryo_width,
        "embryo_center_x" : embryo_center_x,
        "embryo_center_y" : embryo_center_y
    }
    return bounding_roi

bounding_roi = bounding_roi_and_params_from_embryo_mask(mask)
mask.setRoi(bounding_roi["bounding_roi_rect"])
print(bounding_roi)
# File(label='Choose a test image', style='file') image
#  File(label='crop template', style='file') crop_template
#@ DatasetService ds
#@ DatasetIOService io
#@ OpService ops
#@ UIService ui
#@ ImageJ ij
#@ ConvertService convert


from distutils.dir_util import mkpath
import imp
import math
import os
import re
import fnmatch
import json
import logging
from datetime import datetime

from java.io import File
from ij.io import FileSaver
from ij.io import OpenDialog
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin.filter import RankFilters
from fiji.threshold import Auto_Threshold
from ij.plugin.filter import ParticleAnalyzer
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from fiji.selection import Select_Bounding_Box
from ij.plugin import RoiRotator
from ij.plugin import ZProjector
from ij.io import Opener
from ij.plugin import RoiReader
from ij.gui import PointRoi, RotatedRectRoi



def get_rotated_rect_roi_width():
	img = WindowManager.getCurrentImage()
	roi = img.getRoi()
	print roi.getClass()
	x1 = roi.getPolygon().xpoints[0]
	x2 = roi.getPolygon().xpoints[1]
	y1 = roi.getPolygon().ypoints[0]
	y2 = roi.getPolygon().ypoints[1]
	width = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	width = (width / 10 ) * 10 
	print(x1, y1, x2, y2)
	print width


def get_polygon_roi_angle(roi):
	"""Returns an angle between first line in the PolygonRoi and horizontal line

	Args:
		roi (PolygonRoi): intended for rectangular Rois

	Returns:
		double: angle between first line in the PolygonRoi and horizontal line
	"""
	x1 = roi.getPolygon().xpoints[0]
	x2 = roi.getPolygon().xpoints[1]
	y1 = roi.getPolygon().ypoints[0]
	y2 = roi.getPolygon().ypoints[1]
	angle = roi.getAngle(x1, y1, x2, y2)
	# TODO: Find a better way to find angle of rotation of the box
	if angle > 135:
		angle = angle - 180
	if angle < -135:
		angle = 180 - angle
	# logging.info("\tBounding box x-coord:%s, y-coord:%s, rot-angle:%s" %
	# 			 (roi.getPolygon().xpoints, roi.getPolygon().ypoints, angle))
	return angle


def get_rotated_rect_roi_dims(roi):
	x1,	y1,	x2,	y2,	width = roi.getParams()
	height = math.floor(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2))
	# logging.info("Determined rectangular roi height before rounding: %s" % height)
	# height = round(height / 10) * 10
	# logging.info("\tRectangular roi height: %s" % height)
	print("H, W:", (height, width))
	return (width, height)


def midpoint(x1, y1, x2, y2):
    return ((x1 + x2) / 2, (y1 + y2) / 2)

def polygon_to_rotated_rect_roi(roi):
	if isinstance(roi, RotatedRectRoi):
		return roi
	x1 = roi.getPolygon().xpoints[0]
	x2 = roi.getPolygon().xpoints[1]
	x3 = roi.getPolygon().xpoints[2]
	x4 = roi.getPolygon().xpoints[3]
	y1 = roi.getPolygon().ypoints[0]
	y2 = roi.getPolygon().ypoints[1]
	y3 = roi.getPolygon().ypoints[2]
	y4 = roi.getPolygon().ypoints[3]
	dist1 = math.floor(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2))
	dist2 = math.floor(math.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2))
	if dist1 > dist2:
		width = dist2
		rx1, ry1 = midpoint(x1, y1, x4, y4)
		rx2, ry2 = midpoint(x2, y2, x3, y3)
		rot_rect_roi = RotatedRectRoi(rx1, ry1, rx2, ry2, width)
		print(x1, y1, x2, y2, width)
	else:
		width = dist1
		rx1, ry1 = midpoint(x1, y1, x2, y2)
		rx2, ry2 = midpoint(x3, y3, x4, y4)
		rot_rect_roi = RotatedRectRoi(rx1, ry1, rx2, ry2, width)
		print(x2, y2, x3, y3, width)
	return rot_rect_roi

def check_manual_crop_box():
	img = WindowManager.getCurrentImage()
	roi = img.getRoi()
	roi = polygon_to_rotated_rect_roi(roi)
	img.setRoi(roi)
	cropped = img.crop()
	IJ.run(cropped, "Select All", "")

	IJ.run(cropped, "Rotate... ",
            "angle=%s grid=1 interpolation=Bilinear" % round(get_polygon_roi_angle(roi), ndigits=1))
	final_center_x = cropped.getWidth() / 2
	final_center_y = cropped.getHeight() / 2
	box_width, box_height = get_rotated_rect_roi_dims(roi)
	IJ.run(cropped, "Specify...", "width=%s height=%s x=%s y=%s centered" %
	       (box_height, box_width, final_center_x, final_center_y))
	cropped = cropped.crop()
	IJ.run(cropped, "Rotate 90 Degrees Left", "")
	cropped.show()


# get_rotated_rect_roi_width()
check_manual_crop_box()	
	



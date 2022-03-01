
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
	# print roi.getClass()
	x1 = roi.getPolygon().xpoints[0]
	x2 = roi.getPolygon().xpoints[1]
	y1 = roi.getPolygon().ypoints[0]
	y2 = roi.getPolygon().ypoints[1]
	width = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	width = (width / 10) * 10
	# print(x1, y1, x2, y2)
	# print width


def get_polygon_roi_angle(roi):
	"""Returns an angle between longer line in the rectangular PolygonRoi and horizontal line

	Args:
		roi (PolygonRoi): intended for rectangular Rois

	Returns:
		double: angle between longer line in the PolygonRoi and horizontal line
	"""
	x1 = roi.getPolygon().xpoints[0]
	x2 = roi.getPolygon().xpoints[1]
	x3 = roi.getPolygon().xpoints[2]
	y1 = roi.getPolygon().ypoints[0]
	y2 = roi.getPolygon().ypoints[1]
	y3 = roi.getPolygon().ypoints[2]
	dist1 = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	dist2 = math.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)
	if dist1 > dist2:
		if y1 > y2:
			angle = roi.getAngle(x1, y1, x2, y2)
		else:
			angle = roi.getAngle(x2, y2, x1, y1)
	else:
		if y2 > y3:
			angle = roi.getAngle(x2, y2, x3, y3)
		else:
			angle = roi.getAngle(x3, y3, x2, y2)
	if angle > 90:
		angle -= 180
	# logging.info("\tBounding box x-coord:%s, y-coord:%s, rot-angle:%s" %
	# 			 (roi.getPolygon().xpoints, roi.getPolygon().ypoints, angle))
	return angle


def get_rotated_rect_roi_dims(roi):
	x1,	y1,	x2,	y2,	width = roi.getParams()
	height = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	# logging.info("Determined rectangular roi height before rounding: %s" % height)
	# height = round(height / 10) * 10
	# logging.info("\tRectangular roi height: %s" % height)
	# print("H, W:", (height, width))
	return (int(round(width)), int(round(height)))


def midpoint(x1, y1, x2, y2):
	x1, y1, x2, y2 = float(x1), float(y1), float(x2), float(y2)
	return ((x1 + x2) / 2, (y1 + y2) / 2)


def polygon_to_rotated_rect_roi(roi):
	if isinstance(roi, RotatedRectRoi):
		return roi
	print roi.getFloatPolygon().xpoints
	print roi.getFloatPolygon().ypoints
	x1 = roi.getFloatPolygon().xpoints[0]
	x2 = roi.getFloatPolygon().xpoints[1]
	x3 = roi.getFloatPolygon().xpoints[2]
	x4 = roi.getFloatPolygon().xpoints[3]
	y1 = roi.getFloatPolygon().ypoints[0]
	y2 = roi.getFloatPolygon().ypoints[1]
	y3 = roi.getFloatPolygon().ypoints[2]
	y4 = roi.getFloatPolygon().ypoints[3]
	if roi.getNCoordinates() > 4:
		x5 = roi.getFloatPolygon().xpoints[4]
		y5 = roi.getFloatPolygon().ypoints[4]
		if abs(x5 - x1) + abs(y5 - y1) > 1e-05:
			# For some reason after several rotations rectangle polygon can have 5 points
			raise Exception(
				"polygon_to_rotated_rect_roi Can only convert rectangles. Received polygon with %s points. And first and last points are not the same." % roi.getNCoordinates())
	# print("Original polygon: ", x1, y1, x2, y2, x3, y3, x4, y4)
	dist1 = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	dist2 = math.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)
	if dist1 > dist2:
		width = dist2
		rx1, ry1 = midpoint(x1, y1, x4, y4)
		rx2, ry2 = midpoint(x2, y2, x3, y3)
		rot_rect_roi = RotatedRectRoi(rx1, ry1, rx2, ry2, width)
		# print(x1, y1, x2, y2, width)
	else:
		width = dist1
		rx1, ry1 = midpoint(x1, y1, x2, y2)
		rx2, ry2 = midpoint(x3, y3, x4, y4)
		rot_rect_roi = RotatedRectRoi(rx1, ry1, rx2, ry2, width)
		# print(x2, y2, x3, y3, width)
	return rot_rect_roi


def check_manual_crop_box():
	img = WindowManager.getCurrentImage()
	roi_input = img.getRoi()

	# Specify angle of rotation here
	roi_input = RoiRotator.rotate(roi_input, 0)

	roi = polygon_to_rotated_rect_roi(roi_input)
	img.setRoi(roi)
	#print("Dims: ", get_rotated_rect_roi_dims(roi))
	#print("Angle: ", get_polygon_roi_angle(roi))
	# exit()
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

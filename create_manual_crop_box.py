#@ Float   (label="How much to rotate the bounding box? (degrees)", value=0, persist=false, style="slider,format:0.0", min=-90, max=90, stepSize=1) rotation_angle
#@ String (label="Which direction to rotate the embryo", choices={"Right", "Left"}, value=0, style="radioButtonHorizontal") rotation_direction

# Written by Artemiy Golden on Jan 2022 at AK Stelzer Group at Goethe Universitaet Frankfurt am Main
# For detailed documentation go to https://github.com/artgolden/fiji_scripts

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
	x1 = roi.getPolygon().xpoints[0]
	x2 = roi.getPolygon().xpoints[1]
	y1 = roi.getPolygon().ypoints[0]
	y2 = roi.getPolygon().ypoints[1]
	width = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	width = (width / 10) * 10



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
	return angle


def get_rotated_rect_roi_dims(roi):
	x1,	y1,	x2,	y2,	width = roi.getParams()
	height = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	return (int(round(width)), int(round(height)))


def midpoint(x1, y1, x2, y2):
	x1, y1, x2, y2 = float(x1), float(y1), float(x2), float(y2)
	return ((x1 + x2) / 2, (y1 + y2) / 2)


def polygon_to_rotated_rect_roi(roi):
	if isinstance(roi, RotatedRectRoi):
		return roi
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
	dist1 = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	dist2 = math.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)
	if dist1 > dist2:
		width = dist2
		rx1, ry1 = midpoint(x1, y1, x4, y4)
		rx2, ry2 = midpoint(x2, y2, x3, y3)
		rot_rect_roi = RotatedRectRoi(rx1, ry1, rx2, ry2, width)
	else:
		width = dist1
		rx1, ry1 = midpoint(x1, y1, x2, y2)
		rx2, ry2 = midpoint(x3, y3, x4, y4)
		rot_rect_roi = RotatedRectRoi(rx1, ry1, rx2, ry2, width)
	return rot_rect_roi


def check_manual_crop_box():
	img = WindowManager.getCurrentImage()
	roi_input = img.getRoi()

	# Specify angle of rotation here
	roi_input = RoiRotator.rotate(roi_input, round(rotation_angle, ndigits=1))

	roi = polygon_to_rotated_rect_roi(roi_input)
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
	IJ.run(cropped, "Rotate 90 Degrees %s" % rotation_direction, "")
	cropped.show()


check_manual_crop_box()

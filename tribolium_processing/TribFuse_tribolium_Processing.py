#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert
# @ File(label='Choose a directory with datasets', style='directory') datasets_dir
# @ File(label='Choose a file with metadata about embryo directions', style='file') metadata_file
# @ String(label='Dataset prefix', value='MGolden2022A-') dataset_name_prefix
# @ Boolean (label='Compress images?', value=true) compress_on_save
# @ Boolean (label='B-branch: Use previously cropped stacks (if present)?', value=false) use_cropped_cache
# @ Boolean (label='Fusion-branch: Use previously created dataset (works only if everything up to fusion start is done)?', value=false) use_fusion_dataset_cache
# @ Boolean (label='Do B-branch processing?', value=false) do_b_branch
# @ Boolean (label='Do Fusion-branch processing?', value=false) do_fusion_branch
# @ Boolean (label='Do histogram matching adjustment?', value=true) do_histogram_matching
# @ Integer (label='Percentage of overexposed pixels during histogram contrast adjustment', value=1) PERCENT_OVEREXPOSED_PIXELS
# @ Boolean (label='Run fusion-branch only up to downsampled preview', value=true) run_fusion_up_to_downsampled_preview
# @ Boolean (label='Run fusion-branch only up to start of deconvolution', value=False) run_fusion_up_to_start_of_deconvolution
# @ String (choices={"Only Fuse full dataset (BigStitcher)", "Content-based fusion (BigStitcher)", "Fuse+Deconvolve (BigStitcher)", "Deconvolve->Fuse content-based (CLIJx GPU)"}, style="listBox", value="Fuse+Deconvolve") Only_fuse_or_deconvolve
# @ File(label='Choose a cache directory (preferably on an SSD or Ramdisk, otherwise you get no benefits.) (optional param.)', style='directory', required=false) caching_dir
# @ String (label='Thresholding for embryo segmentation', choices={"triangle", "minerror", "mean", "huang2"}, style="radioButtonHorizontal", value="triangle") segmentation_threshold_type
# @ Float (label='Pixel distance X axis for calibration in um', value=1.0, style="format:#.00") pixel_distance_x
# @ Float (label='Pixel distance Y axis for calibration in um', value=1.0, style="format:#.00") pixel_distance_y
# @ Float (label='Pixel distance Z axis for calibration in um', value=4.0, style="format:#.00") pixel_distance_z
# @ Integer (label='Number of deconvolution iterations', value=5) num_deconv_iterations
#  Integer (label='How many times downsample images in XY for preview/segmentation fusion', value=4) image_downsampling_factor

# Written by Artemiy Golden on Jan 2022 at AK Stelzer Group at Goethe Universitaet Frankfurt am Main
# For detailed documentation go to https://github.com/artgolden/tribolium_processing or read the README.md
# Last manual update of this line 2022.5.18 :)

from collections import namedtuple
from copy import deepcopy
from distutils.dir_util import mkpath
import math
import os
import platform
import fnmatch
import json
import logging
from datetime import datetime
import subprocess
import time
import traceback
import sys
import uuid
import timeit




from java.lang.System import getProperty
from java.lang import Runtime
from java.io import File

from ij.io import FileSaver
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin.filter import RankFilters
from fiji.threshold import Auto_Threshold
from ij.plugin.filter import ParticleAnalyzer
from ij.measure import Measurements
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from ij.plugin import RoiRotator
from ij.plugin import Slicer
from ij.plugin import StackCombiner
from ij.io import RoiDecoder

sys.path.append(os.path.join(getProperty("fiji.dir"), "plugins", "tribolium_processing")) 
from tribolium_bigstitcher_utils import *
from tribolium_image_utils import *
from tribolium_math_utils import *
from tribolium_file_utils import *
import tribolium_image_utils

# A very clunky way of importing services for image_utils module. I could find a way to import them directly there, so I have to set the ones that I get from ImageJ2 Script Parameters
tribolium_image_utils.opsServiceImageUtilsLocal = ops
tribolium_image_utils.dsServiceImageUtilsLocal = ds
tribolium_image_utils.convertServiceImageUtilsLocal = convert

# TODO: 
# - Write current limitations and implicit assumptions in the documentation
# - test 2 channel case
# - make option to redo the downsampled images for preview, or better, check if they are the same dimension as planned
# - make an option to load the PSF from a file
# - Since we have to use same min when creating the full dataset, bigstitcher openes every timepoints and checks min/max of the images to determine global one. 
# this is very slow. Probably need to switch to some other way of manually defining min/max.
# - Write install documentation
# - Handle exception when embryo could not be segmented Fusion-branch
# - Generate max time projections uncropped before trying to find the cropbox, so you have projection to define the manual cropbox on
# - Check whether B-branch uses max projections stack cache
# - Convert timepoints for fusion range specification to "free style" string in JSON
# - Handle B-branch exception when embryo could not be segmented:
	# [ERROR] Traceback (most recent call last):
	#   File "tribolium_processing/Tribolium_Processing.py", line 1741, in <module>
	# java.lang.IllegalArgumentException: java.lang.IllegalArgumentException: Row out of range
# - Do more proper checking of the already made PSFs for fusion, need to check that all have been assigned
# - Add a lot more logging to fusion branch. log all bigstitcher commands?
# - restart and log stdout from CLIJ fusion process if it has failed 1 time
# - make Process spawning work for Linux as well
# - Add checking that all raw aligned stacks have been moved from cache dir


EXAMPLE_JSON_METADATA_FILE = """
Example JSON metadata file contents:
{
	"datasets": [
		{
			"ID": 1,
			"channel_1": {
				"specimens_for_directions_1234": [
				5,
				6,
				4,
				7
				]				
			},
			"head_direction": "right",
			"use_manual_bounding_box": false,
			"timepoints_range_to_fuse": [1,100]
		},
		{
			"ID": 3,
			"channel_1": {
				"specimens_for_directions_1234": [
					5,
					6,
					4,
					7
				]
			},
			"channel_2": {
				"specimens_for_directions_1234": [
					0,
					2,
					1,
					3
				]
			},
			"head_direction": "left",
			"use_manual_bounding_box": true,
			"planes_to_keep_per_direction": [
			{
				"start": 1,
				"end": 150
			},
			{
				"start": 10,
				"end": 160
			},
			{
				"start": 1,
				"end": 150
			},
			{
				"start": 1,
				"end": 150
			}
		]
		}
	]
}"""

METADATA_DIR_NAME = "(B1)-Metadata"
TSTACKS_DIR_NAME = "(B3)-TStacks-ZM"
RAW_IMAGES_DIR_NAME = "(P0)-ZStacks-Raw"
RAW_ALIGNED_DIR_NAME = "(P1)-ZStack-AlignedCropped"
RAW_CROPPED_DIR_NAME = "(B2)-ZStacks"
CONTRAST_DIR_NAME = "(B4)-TStacks-ZN"
MONTAGE_DIR_NAME = "(B5)-TStacks-ZN-Montage"
LINK_DIR_FOR_FUSION = "(FUSION)-linked-dataset"

B_BRANCH_ERROR_FILE_NAME = "B_BRANCH_ERRORED"
B_BRANCH_FINISHED_FILE_NAME = "B_BRANCH_FINISHED"
B_BRANCH_LOCK_FILE_NAME = "B_BRANCH_LOCK"

FUSION_BRANCH_ERROR_FILE_NAME = "FUSION_BRANCH_ERRORED"
FUSION_BRANCH_FINISHED_FILE_NAME = "FUSION_BRANCH_FINISHED"
FUSION_BRANCH_LOCK_FILE_NAME = "FUSION_BRANCH_LOCK"

MANUAL_CROP_BOX_FILE_NAME = "manual_crop_box"


DEFAULT_CROP_BOX_WIDTH = 1100
MINIMUM_CROP_BOX_WIDTH = 1000


# Percentage of overexposed pixels during histogram contrast adjustment
# PERCENT_OVEREXPOSED_PIXELS = 1

opsServiceImageUtilsLocal = opsServiceImageUtilsLocal
dsServiceImageUtilsLocal = dsServiceImageUtilsLocal
convertServiceImageUtilsLocal = convertServiceImageUtilsLocal

# Reassigning Scripting parameters to get warnings about missing variables only here

ops = ops
ds = ds
convert = convert

USE_CROPPED_CACHE = use_cropped_cache
USE_FUSION_DATASET_CACHE = use_fusion_dataset_cache
DO_HISTOGRAM_MATCHING = do_histogram_matching
PERCENT_OVEREXPOSED_PIXELS = PERCENT_OVEREXPOSED_PIXELS
COMPRESS_ON_SAVE = compress_on_save
DO_B_BRANCH = do_b_branch
DO_FUSION_BRANCH = do_fusion_branch
RUN_FUSION_UP_TO_DOWNSAMPLED_PREVIEW = run_fusion_up_to_downsampled_preview
RUN_FUSION_BRANCH_UP_TO_START_OF_FUSION = run_fusion_up_to_start_of_deconvolution
CACHING_DIR = os.path.join(datasets_dir.getAbsolutePath(), "tmp")
if not os.path.exists(CACHING_DIR):
	mkpath(CACHING_DIR)
if not caching_dir is None: 
	CACHING_DIR = caching_dir.getAbsolutePath()
IMAGE_DOWNSAMPLING_FACTOR_XY = image_downsampling_factor = 4
ONLY_FUSE_FULL = False
FUSE_AND_DECONVOLVE_FULL = False
FUSE_CONTENT_BASED = False
FUSE_CLIJ_GPU = False
SEGMENTATION_THRESHOLD_TYPE = segmentation_threshold_type

IMAGE_PIXEL_DISTANCE_X = pixel_distance_x
IMAGE_PIXEL_DISTANCE_Y = pixel_distance_y
IMAGE_PIXEL_DISTANCE_Z = pixel_distance_z

NUM_DECONV_ITERATIONS = num_deconv_iterations

FUSED_DIR_NAME = "(F0)-GenericFusion"
if Only_fuse_or_deconvolve == "Only Fuse full dataset (BigStitcher)":
	ONLY_FUSE_FULL = True
	FUSED_DIR_NAME = "(F0)-ZStacks-SimpleFusion-ZS"
if Only_fuse_or_deconvolve == "Content-based fusion (BigStitcher)":
	FUSE_CONTENT_BASED = True
	FUSED_DIR_NAME = "(F0)-ZStacks-ContentBasedFusion-ZS"
if Only_fuse_or_deconvolve == "Fuse+Deconvolve (BigStitcher)":
	FUSE_AND_DECONVOLVE_FULL = True
	FUSED_DIR_NAME = "(F0)-ZStacks-BigStitcherDeconFusion-ZS"
if Only_fuse_or_deconvolve == "Deconvolve->Fuse content-based (CLIJx GPU)":
	FUSE_CLIJ_GPU = True
	FUSED_DIR_NAME = "(F0)-ZStacks-DeconContentBasedFusion-ZS"

CANVAS_SIZE_FOR_EMBRYO_PREVIEW = 1400 / IMAGE_DOWNSAMPLING_FACTOR_XY

DEBUG_USE_FUSED_PREVIEW_CACHE = False

class DatasetMetadata:
	id = 0
	number_of_directions = 0
	number_of_channels = 0
	raw_images_dir_path = ""
	root_dir = ""
	datasets_dir = ""
	specimen_directions_in_channels = None
	embryo_head_direction = ""
	use_manual_b_branch_bounding_box = False
	planes_to_keep_per_direction_b_branch = None
	tp_fuse_start = 1
	tp_fuse_end = -1
	raw_aligned_stacks_dir = ""

	def __init__(
				self, 
				id, 
				root_dir, 
				specimen_directions_in_channels, 
				embryo_head_direction,
				use_manual_b_branch_bounding_box=False,
				planes_to_keep_per_direction_b_branch=None,
				tp_fuse_start=1,
				tp_fuse_end=-1):
		self.id = id
		self.root_dir = root_dir
		self.datasets_dir = os.path.dirname(root_dir)
		self.embryo_head_direction = embryo_head_direction
		self.specimen_directions_in_channels = specimen_directions_in_channels
		self.use_manual_b_branch_bounding_box = use_manual_b_branch_bounding_box
		self.planes_to_keep_per_direction_b_branch = planes_to_keep_per_direction_b_branch
		self.number_of_channels = len(specimen_directions_in_channels)
		self.number_of_directions = len(specimen_directions_in_channels[0])
		self.raw_images_dir_path = os.path.join(root_dir, RAW_IMAGES_DIR_NAME)
		self.tp_fuse_start = tp_fuse_start
		self.tp_fuse_end = tp_fuse_end
		self.raw_aligned_stacks_dir = os.path.join(self.root_dir, RAW_ALIGNED_DIR_NAME)


def get_DatasetMetadata_obj_from_metadata_dict(metadata_dict, datasets_dir):

	dataset_root_dir = os.path.dirname(get_raw_images_dir(datasets_dir, metadata_dict["ID"]))

	specimen_directions_in_channels = []
	for chan in range(3):
		channel = "channel_%s" % chan
		if channel in metadata_dict:
			specimen_directions_in_channels.append(
				tuple(metadata_dict[channel]["specimens_for_directions_1234"]))
	
	specimen_directions_in_channels = tuple(specimen_directions_in_channels)

	if "planes_to_keep_per_direction" in metadata_dict:
		planes_to_keep_per_direction_b_branch=metadata_dict["planes_to_keep_per_direction"]
	else:
		planes_to_keep_per_direction_b_branch = None
	
	tp_fuse_start = 1
	tp_fuse_end = -1
	if "timepoints_range_to_fuse" in metadata_dict.keys():
		tp_fuse_start = metadata_dict["timepoints_range_to_fuse"][0]
		tp_fuse_end = metadata_dict["timepoints_range_to_fuse"][1]

	return DatasetMetadata(
							id=metadata_dict["ID"], 
							root_dir=dataset_root_dir, 
							specimen_directions_in_channels=specimen_directions_in_channels, 
							embryo_head_direction = metadata_dict["head_direction"],
							use_manual_b_branch_bounding_box=metadata_dict["use_manual_bounding_box"],
							planes_to_keep_per_direction_b_branch=planes_to_keep_per_direction_b_branch,
							tp_fuse_start=tp_fuse_start,
							tp_fuse_end=tp_fuse_end,
							)
		
SegmentationResults = namedtuple("SegmentationResults", ["transformation_embryo_to_center", "rotation_angle_z", "rotation_angle_y", "embryo_crop_box"])


# 						Generic functions

def parse_range_string(astr):
	result = set()
	for part in astr.split(','):
		x = part.split('-')
		result.update(range(int(x[0]), int(x[-1]) + 1))
	return sorted(result)

def get_free_memory_in_GB():
	free_mem = round(float(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1024 / 1024 / 1024, 3)
	return free_mem

def start_moving_files(file_paths, destination_dir):
	for path in file_paths:
		logging.info("Started move process of the file %s to %s" % (path, destination_dir))
		move_command = ["move", path, destination_dir]
		p = subprocess.Popen(move_command, shell=True)


#						B-branch functions

def create_crop_template(max_time_projection, meta_dir, dataset_metadata_obj, dataset_minimal_crop_box_dims, use_dataset_box_dims):
	"""Crop max_time_projection, create crop template .roi object, save it in meta_dir.

	Args:
		max_time_projection (ImagePlus): a single image of max projection of a time stack of max projections
		meta_dir (str): full path to metadata dir
		dataset (DatasetMetadata): metadata about the dataset to extract information about embryo orientation
		dataset_maximal_crop_box_dims (int, int): width, heigth
		use_dataset_box_dims (bool):

	Returns:
		Roi: bounding box around the embryo with proper rotation and size
		ImagePlus: cropped max time projection
		(int, int): updated dataset maximal crop box dims
	"""
	updated_dataset_maximal_crop_box_dims = dataset_minimal_crop_box_dims
	imp = max_time_projection

	if dataset_metadata_obj.use_manual_b_branch_bounding_box == True:
		logging.info("\tUsing manualy specified crop box.")
		manual_roi = RoiDecoder(os.path.join(meta_dir, "%s.roi" % MANUAL_CROP_BOX_FILE_NAME)).getRoi()
		rot_angle = round(get_rectangular_polygon_roi_angle(manual_roi), ndigits=1)
		# Always assuming that embryo is horizontal
		embryo_length, embryo_width = get_rotated_rect_polygon_roi_dims(manual_roi)
		embryo_center_x = manual_roi.getContourCentroid()[0]
		embryo_center_y = manual_roi.getContourCentroid()[1]
	else:
		ip = imp.getProcessor().duplicate()  # get pixel array?, as a copy
		mask = ImagePlus("mask", ip)
		radius = 4
		RankFilters().rank(ip, radius, RankFilters.MEDIAN)
		hist = ip.getHistogram()
		triag_threshold = Auto_Threshold.Triangle(hist)
		ip.setThreshold(triag_threshold, float("inf"), ImageProcessor.NO_LUT_UPDATE)
		IJ.run(mask, "Convert to Mask", "")
		# Repeated Erosion and Dilation to remove dirt, which shows as protrusions on the embryo mask
		# Using IJ.run() for erosion and dilation has some persistent state after run
		# that produces inconsistent behaviour and bugs
		# This is why I am using ByteProcessor methods directly
		smoothen_embryo_mask(mask, 40)
		IJ.run("Clear Results")
		table = ResultsTable()
		# Initialise without display (boolean value is ignored)
		roim = RoiManager(False)
		MIN_PARTICLE_SIZE = 10000  # pixel ^ 2
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
		logging.info("\t Angle of the elipse fitted onto the embryo: %s", rot_angle)
		roi_arr = roim.getRoisAsArray()
		roim.runCommand('reset')
		imp.setRoi(roi_arr[len(roi_arr) - 1])
		IJ.run(imp, "Select Bounding Box (guess background color)", "")
		bounding_roi = imp.getRoi()
		# Crop the image so the borders are equidistant from the borders of the bounding box
		# for the center of rotation to be aligned roughly with the center of the embryo
		bounding_rectangle = bounding_roi.getBounds()
		# extracting width field from java.awt.Rectangle
		i = bounding_rectangle.x
		j = bounding_rectangle.y
		We = bounding_rectangle.width
		He = bounding_rectangle.height
		W = imp.getWidth()
		H = imp.getHeight()
		ax = i
		bx = W - i - We
		ay = j
		by = H - j - He
		xn = i - min(ax, bx)
		yn = j - min(ay, by)
		Wn = We + 2 * min(ax, bx)
		Hn = He + 2 * min(ay, by)
		IJ.run(imp, "Specify...", "width=%s height=%s x=%s y=%s" % (Wn, Hn, xn, yn))
		equidistant_crop_roi = imp.getRoi()
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

		# Extract the precise embryo center in the rotated cropped image
		roim = RoiManager(False)
		MIN_PARTICLE_SIZE = 10000  # pixel ^ 2
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
		IJ.run(mask_cropped, "Select Bounding Box (guess background color)", "")
		cropped_bounding_roi = mask_cropped.getRoi()
		embryo_length = cropped_bounding_roi.getBounds().width
		embryo_width = cropped_bounding_roi.getBounds().height
		logging.info("\tEmbryo dims: (%s, %s)" % (embryo_length, embryo_width))
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
	
	logging.info("\tEmbryo center: (%s, %s)" % (embryo_center_x, embryo_center_y))

	if use_dataset_box_dims == True:
		dataset_width, dataset_height = dataset_minimal_crop_box_dims
		logging.info("\tUsing dataset's global value of the crop box dims: %s, %s to create the crop template" % (dataset_width, dataset_height))
		box_width, box_height = dataset_width, dataset_height
		IJ.run(imp, "Specify...", "width=%s height=%s x=%s y=%s centered" %
					(dataset_width, dataset_height, embryo_center_x, embryo_center_y))
	else:
		dataset_width, dataset_height = dataset_minimal_crop_box_dims
		box_width = embryo_length + 12 + 4 - embryo_length % 4
		box_height = embryo_width + 12 + 4 - embryo_width % 4
		box_width = max(box_width, dataset_width)
		box_height = max(box_height, dataset_height)
		IJ.run(imp, "Specify...", "width=%s height=%s x=%s y=%s centered" %
					(box_width, box_height, embryo_center_x, embryo_center_y))
		updated_dataset_maximal_crop_box_dims = (box_width, box_height)
	bounding_roi = imp.getRoi()
	bounding_roi_rot = RoiRotator.rotate(bounding_roi, -rot_angle)

	# Check if the rotated bounding box extend over the image
	if is_polygon_roi_overlapping_image_edges(imp, bounding_roi_rot) == True:
		if use_dataset_box_dims == True:
			raise Exception(
				"Fixed bounding box width for the whole dataset cannot be applied to all dimensions.")
		else:
			raise Exception(
				"Could not find a bounding box that contains whole embryo and does not protrude out of the image.")
	if updated_dataset_maximal_crop_box_dims != dataset_minimal_crop_box_dims:
		logging.info(
					"\tUpdated the global dataset min crop box dims from %s to: %s" % (dataset_minimal_crop_box_dims, updated_dataset_maximal_crop_box_dims))

	imp.setRoi(bounding_roi_rot, True)
	crop_template = imp.getRoi()
	if dataset_metadata_obj.use_manual_b_branch_bounding_box == False:
		roim.runCommand('reset')
		roim.addRoi(crop_template)
		roim.select(0)
		roim.runCommand("Save", os.path.join(meta_dir, "crop_template.roi"))
		roim.close()
	# This crops by a unrotated bounding box created around the rotated selection box
	final_imp = imp.crop() # so later, when we rotate and crop again there will be no blacked corners in the image.
	IJ.run(final_imp, "Select All", "")
	IJ.run(final_imp, "Rotate... ", "angle=%s grid=1 interpolation=Bilinear" % rot_angle)
	final_center_x = final_imp.getWidth() / 2
	final_center_y = final_imp.getHeight() / 2
	IJ.run(final_imp, "Specify...", "width=%s height=%s x=%s y=%s centered" %
			(box_width, box_height, final_center_x, final_center_y))
	cropped_max_time_proj = final_imp.crop()
	if dataset_metadata_obj.embryo_head_direction == "left":
		IJ.run(cropped_max_time_proj, "Rotate 90 Degrees Right", "")

		# So that illumination is comming from the right.
		# Right now everything is based on that illumination on the images from mDSLM is always coming from the top.
		IJ.run(cropped_max_time_proj, "Flip Horizontally", "")
	else:
		IJ.run(cropped_max_time_proj, "Rotate 90 Degrees Left", "")
	

	return (crop_template, cropped_max_time_proj, updated_dataset_maximal_crop_box_dims)

def find_planes_to_keep(zstack, meta_dir, manual_planes_to_keep):
	"""Calculates planes to keep, so that the embryo is centered in them. 
	Saves a cropped Y_projected stack for user assessment in the metadata directory. 

	Args:
		zstack (ImagePlus): cropped Z-stack
		meta_dir (str): absolute path to metadata directory to save the Y-max projection of the image for the user to asses the plane selection 

	Returns:
		(int, int): (start_plane, stop_plane) ends have to be included.
	"""
	IJ.run(zstack, "Properties...",
			"pixel_width=1.0000 pixel_height=1.0000 voxel_depth=4.0000")
	for_user_assessment = zstack.duplicate()

	if manual_planes_to_keep is None:
		# This variant of implementation does not recalculate the pixel values and does not expand the image.
		# So later, when we convert to mask, the image will be number_of_planes in height, and not 4*number_of_planes
		# as would be with the Reslice made from GUI.
		# This does not account for voxel_depth, so the resliced image is not rectangular.
		resliced = Slicer.reslice(Slicer(), zstack)

		max_proj_reslice = z_project_a_stack(resliced)
		ip = max_proj_reslice.getProcessor()
		radius = 4
		RankFilters().rank(ip, radius, RankFilters.MEDIAN)

		hist = ip.getHistogram()
		triag_threshold = Auto_Threshold.Mean(hist)
		ip.setThreshold(triag_threshold, float("inf"), ImageProcessor.NO_LUT_UPDATE)
		IJ.run(max_proj_reslice, "Convert to Mask", "")

		table = ResultsTable()
		# Initialise RoiManager without display (boolean value is ignored)
		roim = RoiManager(False)
		MIN_PARTICLE_SIZE = 2000  # pixel ^ 2
		MAX_PARTICLE_SIZE = float("inf")
		ParticleAnalyzer.setRoiManager(roim)
		pa = ParticleAnalyzer(
				ParticleAnalyzer.ADD_TO_MANAGER,
				Measurements.RECT,
				table,
				MIN_PARTICLE_SIZE,
				MAX_PARTICLE_SIZE)
		pa.analyze(max_proj_reslice)
		box_start = table.getValue("BY", 0)
		box_width = table.getValue("Height", 0)
		middle_y = box_start + box_width / 2
		start_plane = int(round(middle_y - 75))
		end_plane = int(round(middle_y + 74))

		if start_plane < 1:
			if start_plane < -5:
				raise Exception(
					"Embryo is more than 5 planes off center in Z-direction, for the 150 planes cropping.")
			start_plane = 1
			end_plane = 150
		if end_plane > zstack.getNSlices():
			if end_plane > zstack.getNSlices() + 5:
				raise Exception(
					"Embryo is more than 5 planes off center in Z-direction, for the 150 planes cropping.")
			end_plane = zstack.getNSlices()
			start_plane = zstack.getNSlices() - 149
	else:
		start_plane, end_plane = manual_planes_to_keep
		logging.info("\t Using manually specified planes to keep.")
		if start_plane <= 0 or end_plane > zstack.getNSlices():
			raise Exception(
				"Manually specified planes to keep cannot be applied to a stack with %s planes." % zstack.getNSlices())
	
	logging.info("\t selected planes to keep: %s-%s" % (start_plane, end_plane))
	logging.info("\t Cropping Y-projection for user assessment with %s planes" %
				 for_user_assessment.getNSlices())

	# Save cropped Y-projection for user assessment
	for_user_assessment = subset_planes(
		for_user_assessment, (start_plane, end_plane))

	#TODO: Find how to make this headless
	IJ.run(for_user_assessment, "Reslice [/]...", "output=1 start=Top")
	for_user_assessment = WindowManager.getCurrentImage()
	for_user_assessment.hide()

	for_user_assessment = z_project_a_stack(for_user_assessment)
	fs = FileSaver(for_user_assessment)
	fs.saveAsTiff(os.path.join(
		meta_dir, "Y_projected_raw_stack_for_assessment_of_plane_selection.tif"))

	return (start_plane, end_plane)


#						Fusion-branch functions

def get_embryo_mask(image_16bit, threshold_type):
	ip = image_16bit.getProcessor().duplicate()  # get pixel array?, as a copy
	mask = ImagePlus("mask", ip)
	radius = int(round(mask.getWidth() * 0.0084))
	RankFilters().rank(ip, radius, RankFilters.MEDIAN)
	hist = ip.getHistogram()
	if threshold_type == "triangle":
		min_threshold = Auto_Threshold.Triangle(hist)
	if threshold_type == "mean":
		min_threshold = Auto_Threshold.Mean(hist)
	if threshold_type == "huang2":
		min_threshold = Auto_Threshold.Huang2(hist)
	if threshold_type == "minerror":
		min_threshold = Auto_Threshold.MinErrorI(hist)
	ip.setThreshold(min_threshold, float("inf"), ImageProcessor.NO_LUT_UPDATE)
	IJ.run(mask, "Convert to Mask", "")
	return mask

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

def bounding_roi_and_params_from_embryo_mask(mask, spacing_around_embryo=10):
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
	# print("\t Angle of the elipse fitted onto the embryo: %s", rot_angle)

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
	# print("\tEmbryo dims: (%s, %s)" % (embryo_length, embryo_width))
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
	
	
	# print("\tEmbryo center: (%s, %s)" % (embryo_center_x, embryo_center_y))
	box_width = embryo_length + spacing_around_embryo + 4 - embryo_length % 4
	box_height = embryo_width + spacing_around_embryo + 4 - embryo_width % 4
	IJ.run(mask, "Specify...", "width=%s height=%s x=%s y=%s centered" %
			(box_width, box_height, embryo_center_x, embryo_center_y))
	bounding_roi = mask.getRoi()
	bounding_roi_rot = RoiRotator.rotate(bounding_roi, -rot_angle)
	mask.setRoi(bounding_roi_rot)
	bounding_roi = {
		"bounding_roi_rect" : bounding_roi_rot,
		"embryo_length" : embryo_length,
		"embryo_width" : embryo_width,
		"embryo_center_x" : embryo_center_x,
		"embryo_center_y" : embryo_center_y,
		"bounding_rect_angle" : rot_angle
	}
	return bounding_roi

def get_embryo_bounding_rectangle_and_params(max_projection, show_mask=False, fill_zeros=False):
	if fill_zeros:
		fill_zero_pixels_with_median_intensity(max_projection)

	mask = get_embryo_mask(image_16bit=max_projection, threshold_type=SEGMENTATION_THRESHOLD_TYPE)


	if show_mask:
		mask.show() 

	num_mask_smoothening_iterations = int(round(mask.getWidth() * 0.028))
	# removing dirt, fixing edges
	smoothen_embryo_mask(mask, num_mask_smoothening_iterations)

	bounding_roi_stuff = bounding_roi_and_params_from_embryo_mask(mask=mask, spacing_around_embryo=8)
	bounding_roi_stuff["mask"] = mask
	return bounding_roi_stuff

def create_and_fuse_preview_dataset(dataset_dir, file_pattern, angles_from_to, view_stack_dims):
	dataset_xml_name = "dataset.xml"
	dataset_xml = os.path.join(dataset_dir, "dataset.xml")
	sigma = 2.0
	threshold = 0.01
	max_detections_per_view = 1000

	redundancy = 1
	significance = 10
	allowed_error_for_ransac = 2
	number_of_ransac_iterations = "Normal"

	reference_timepoint = 1

	output_fused_path = os.path.join(dataset_dir, "fused.xml")



#  
	IJ.run(
		"Define dataset ...",
		"define_dataset=[Manual Loader (TIFF only, ImageJ Opener)] project_filename=%s multiple_timepoints=[NO (one time-point)] multiple_channels=[NO (one channel)] _____multiple_illumination_directions=[NO (one illumination direction)] multiple_angles=[YES (one file per angle)] multiple_tiles=[NO (one tile)] image_file_directory=%s image_file_pattern=%s acquisition_angles_=%s-%s calibration_type=[Same voxel-size for all views] calibration_definition=[User define voxel-size(s)] imglib2_data_container=[ArrayImg (faster)] pixel_distance_x=1.00000 pixel_distance_y=1.00000 pixel_distance_z=1.00000 pixel_unit=um" % (dataset_xml_name, dataset_dir, file_pattern, angles_from_to[0], angles_from_to[1]))

	IJ.run(
		"Detect Interest Points for Pairwise Registration", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] type_of_interest_point_detection=Difference-of-Gaussian label_interest_points=beads limit_amount_of_detections subpixel_localization=[3-dimensional quadratic fit] interest_point_specification=[Advanced ...] downsample_xy=[Match Z Resolution (less downsampling)] downsample_z=1x use_same_min sigma=%s threshold=%s find_maxima maximum_number=%s type_of_detections_to_use=Brightest compute_on=[CPU (Java)]"
		% (dataset_xml, sigma, threshold, max_detections_per_view))

	IJ.run(
		"Register Dataset based on Interest Points",
		"select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] registration_algorithm=[Fast descriptor-based (rotation invariant)] registration_over_time=[Register timepoints individually] registration_in_between_views=[Only compare overlapping views (according to current transformations)] interest_points=beads fix_views=[Fix first view] map_back_views=[Do not map back (use this if views are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 redundancy=%s significance=%s allowed_error_for_ransac=%s number_of_ransac_iterations=%s"
		% (dataset_xml, redundancy, significance, allowed_error_for_ransac, number_of_ransac_iterations))

	box_coords = {
		"x_min" : 0 + 1,
		"y_min" : 0 + 1,
		"z_min" : 0 + 1,
		"x_max" : view_stack_dims[0] - 1,
		"y_max" : view_stack_dims[1] - 1,
		"z_max" : view_stack_dims[2] - 2,
	}
	define_bounding_box_for_fusion(dataset_xml, box_coords, "max_box")

	logging.info("saving fused dataset to: %s" % output_fused_path)
	fuse_dataset_to_tiff(dataset_xml, "max_box", output_fused_path)

def segment_embryo_and_fuse_again_cropping_around_embryo(raw_dataset_xml_path, fused_xml_path, z_projection, y_projection):
	dataset_dir = os.path.dirname(raw_dataset_xml_path)
	cropped_fused_dir = os.path.join(dataset_dir, "cropped_fusion")
	if not os.path.exists(cropped_fused_dir):
		mkpath(cropped_fused_dir)

	z_bounding_roi = get_embryo_bounding_rectangle_and_params(max_projection=z_projection, show_mask=False)
	y_bounding_roi = get_embryo_bounding_rectangle_and_params(max_projection=y_projection, show_mask=False)

	if abs(z_bounding_roi["embryo_center_x"] - y_bounding_roi["embryo_center_x"]) > 5:
		print("Could not reliably determine embryo X center position from Z and Y projections.")
		logging.info("Could not reliably determine embryo X center position from Z and Y projections.")
		return False

	logging.info("Embryo box params from Z projection:\n%s" % z_bounding_roi)
	logging.info("Embryo box params from Y projection:\n%s " % y_bounding_roi)

	
	save_roi(z_bounding_roi["bounding_roi_rect"], os.path.join(dataset_dir, "z_projection_bounding_rectangle.roi"))
	save_roi(y_bounding_roi["bounding_roi_rect"], os.path.join(dataset_dir, "y_projection_bounding_rectangle.roi"))
	save_tiff(z_bounding_roi["mask"], os.path.join(dataset_dir, "z_projection_segmentation_mask.tiff"))
	save_tiff(y_bounding_roi["mask"], os.path.join(dataset_dir, "y_projection_segmentation_mask.tiff"))


	transformation_embryo_to_center = get_4_4_identity_matrix()
	transformation_embryo_to_center[0][3] = z_bounding_roi["embryo_center_x"]
	transformation_embryo_to_center[1][3] = z_bounding_roi["embryo_center_y"]
	transformation_embryo_to_center[2][3] = y_bounding_roi["embryo_center_y"]
	transf_to_raw = get_fusion_tranformation_from_xml_file(fused_xml_path)
	if transf_to_raw == False:
		print("ERROR: could not extract transformation matrix from fusion XML file!")
		logging.info("ERROR: could not extract transformation matrix from fusion XML file!")
		exit(1)
	transformation_embryo_to_center = multiply_matrices(transformation_embryo_to_center, transf_to_raw)
	transformation_embryo_to_center = inverse_matrix(transformation_embryo_to_center)

	apply_transformation_bigstitcher_dataset_0_timepoint(raw_dataset_xml_path, transformation_embryo_to_center)
	embryo_crop_box = {
		"x_min" : -1 * int(z_bounding_roi["embryo_length"] / 2 + 5),
		"y_min" : -1 * int(z_bounding_roi["embryo_width"] / 2 + 5),
		"z_min" : -1 * int(y_bounding_roi["embryo_width"] / 2 + 5),
		"x_max" : int(z_bounding_roi["embryo_length"] / 2 + 5),
		"y_max" : int(z_bounding_roi["embryo_width"] / 2 + 5),
		"z_max" : int(y_bounding_roi["embryo_width"] / 2 + 5)
	}

	rotation_angle_z = round(z_bounding_roi["bounding_rect_angle"], 1)
	rotation_angle_y = -1 * round(y_bounding_roi["bounding_rect_angle"], 1)
	rotate_bigstitcher_dataset_0_timepoint(raw_dataset_xml_path, "z", rotation_angle_z)
	rotate_bigstitcher_dataset_0_timepoint(raw_dataset_xml_path, "y", rotation_angle_y)


	define_bounding_box_for_fusion(raw_dataset_xml_path, embryo_crop_box, "embryo_cropped")
	cropped_fused_path = os.path.join(cropped_fused_dir, "cropped_fusion.xml")
	cropped_fusion_stack_path = os.path.join(cropped_fused_dir, "fused_tp_0_ch_0.tif")
	fuse_dataset_to_tiff(raw_dataset_xml_path, "embryo_cropped", cropped_fused_path)
	zy_cropped_projections = Z_Y_projection_montage_from_image_stack(IJ.openImage(cropped_fusion_stack_path))

	segmentation_results = SegmentationResults(transformation_embryo_to_center, rotation_angle_z, rotation_angle_y, embryo_crop_box)

	return (zy_cropped_projections, segmentation_results)

def upscale_bounding_box(bounding_box, downsampling_factor_XY, pixel_sizes):
	z_scaling_factor = pixel_sizes[2] / pixel_sizes[0]
	upscaled_box = {
		"x_min" : bounding_box["x_min"] * downsampling_factor_XY,
		"y_min" : bounding_box["y_min"] * downsampling_factor_XY,
		"z_min" : bounding_box["z_min"] * z_scaling_factor,
		"x_max" : bounding_box["x_max"] * downsampling_factor_XY,
		"y_max" : bounding_box["y_max"] * downsampling_factor_XY,
		"z_max" : bounding_box["z_max"] * z_scaling_factor 
	}
	return upscaled_box

def rotate_bounding_box_vertically(box):
	vertical_box = {
		"x_min" : box["y_min"],
		"y_min" : box["x_min"],
		"z_min" : box["z_min"],
		"x_max" : box["y_max"],
		"y_max" : box["x_max"],
		"z_max" : box["z_max"]
	}
	return vertical_box

def upscale_translation_matrix(matrix, downsampling_factor_XY):
	upscaled_matrix = matrix
	upscaled_matrix[0][3] = upscaled_matrix[0][3] * downsampling_factor_XY
	upscaled_matrix[1][3] = upscaled_matrix[1][3] * downsampling_factor_XY
	upscaled_matrix[2][3] = upscaled_matrix[2][3] * downsampling_factor_XY
	return upscaled_matrix

def create_and_register_full_dataset(dataset_dir, dataset_name, file_pattern, timepoints_list, angles_from_to, calibration_pixel_distances, reference_timepoint=1):
	dataset_xml_name = "%s.xml" % dataset_name
	dataset_xml = os.path.join(dataset_dir, "%s.xml" % dataset_name)
	sigma = 2.0
	threshold = 0.01
	max_detections_per_view = 1000

	redundancy = 1
	significance = 10
	allowed_error_for_ransac = 2
	number_of_ransac_iterations = "Normal"



	define_string = "define_dataset=[Manual Loader (TIFF only, ImageJ Opener)] project_filename=%s multiple_timepoints=[YES (one file per time-point)] multiple_channels=[NO (one channel)] _____multiple_illumination_directions=[NO (one illumination direction)] multiple_angles=[YES (one file per angle)] multiple_tiles=[NO (one tile)] image_file_directory=%s image_file_pattern=%s timepoints_=%s acquisition_angles_=%s-%s calibration_type=[Same voxel-size for all views] calibration_definition=[User define voxel-size(s)] imglib2_data_container=[ArrayImg (faster)] pixel_distance_x=%s pixel_distance_y=%s pixel_distance_z=%s pixel_unit=um" % (dataset_xml_name, 
																					dataset_dir, 
																					file_pattern, 
																					timepoints_list, 
																					angles_from_to[0], 
																					angles_from_to[1],
																					calibration_pixel_distances[0],
																					calibration_pixel_distances[1],
																					calibration_pixel_distances[2]
																					)
	logging.info("Defining the full dataset: " + define_string)
	IJ.run(
		"Define dataset ...", define_string)

	# use_same_min parameter is very important. If you do not use same min max intensity, point detection will not work for low intensity datasets
	detect_string = "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] type_of_interest_point_detection=Difference-of-Gaussian label_interest_points=beads limit_amount_of_detections subpixel_localization=[3-dimensional quadratic fit] interest_point_specification=[Advanced ...] downsample_xy=[Match Z Resolution (less downsampling)] downsample_z=1x use_same_min sigma=%s threshold=%s find_maxima maximum_number=%s type_of_detections_to_use=Brightest compute_on=[CPU (Java)]" % (dataset_xml, sigma, threshold, max_detections_per_view)
	logging.info("Detecting interest points of the full dataset: " + detect_string)
	IJ.run("Detect Interest Points for Pairwise Registration", detect_string)

	register_string = "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] registration_algorithm=[Fast descriptor-based (rotation invariant)] registration_over_time=[Register timepoints individually] registration_in_between_views=[Only compare overlapping views (according to current transformations)] interest_points=beads fix_views=[Fix first view] map_back_views=[Do not map back (use this if views are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 redundancy=%s significance=%s allowed_error_for_ransac=%s number_of_ransac_iterations=%s" % (dataset_xml, redundancy, significance, allowed_error_for_ransac, number_of_ransac_iterations)
	logging.info("Registering individual timepoints of the full dataset: " + register_string)
	IJ.run("Register Dataset based on Interest Points",register_string)

	register_to_1_tp_string = "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] registration_algorithm=[Fast descriptor-based (rotation invariant)] registration_over_time=[Match against one reference timepoint (no global optimization)] registration_in_between_views=[Only compare overlapping views (according to current transformations)] interest_points=beads reference=%s consider_each_timepoint_as_rigid_unit transformation=Translation regularize_model model_to_regularize_with=Rigid lamba=0.10 redundancy=1 significance=10 allowed_error_for_ransac=2 number_of_ransac_iterations=Normal interestpoint_grouping=[Group interest points (simply combine all in one virtual view)] interest=5" % (dataset_xml, reference_timepoint)
	logging.info("Registering timepoints to tp %s of the full dataset: %s" % (reference_timepoint, register_to_1_tp_string))
	IJ.run("Register Dataset based on Interest Points", register_to_1_tp_string)

def link_all_files_from_directions_to_folder(directions_dirs_list, symlinked_folder):
	# Currently unused
	for direction_dir in directions_dirs_list:
		for file_path in get_tiffs_in_directory(direction_dir):
			symlink_path = os.path.join(symlinked_folder, os.path.basename(file_path))
			if not os.path.exists(symlink_path):
				if "Linux" in platform.platform():
					logging.info("Detected Linux operating system, using symlinks for fusion folder.")
					os.symlink(file_path, symlink_path)
				elif "Windows" in platform.platform():
					logging.info("Detected Windows operating system, using hard links for fusion folder.")
					os.system('cmd /c "mklink /H %s %s"' % (symlink_path, file_path))
				else:
					print("ERROR: Could not determine Operating System type. Exiting.")
					logging.info("ERROR: Could not determine Operating System type. Exiting.")
					exit(1)
	
def get_timepoint_list_from_folder(dir_path, tp_start, tp_end):
	timepoints = []
	for file in get_tiffs_in_directory(dir_path):
		filename = FredericFile(os.path.basename(file))
		try:
			timepoints.append(int(filename.time_point))
		except:
			print("ERROR: could not extract timepoint from file name. %s" % file)
			logging.info("ERROR: could not extract timepoint from file name. %s" % file)
			return False
	# removing duplicates and sorting
	timepoints = list(dict.fromkeys(timepoints))
	timepoints.sort()

	# subset timepoints
	if tp_end != -1:
		timepoints = [timepoints[i] for i in range(timepoints.index(tp_start), timepoints.index(tp_end) + 1)]  
	return timepoints 

def start_deconv_fuse_timepoint_process(raw_transformed_paths, psf_paths, fused_output_path, num_deconv_iter):
	fiji_exe = os.path.join(getProperty("fiji.dir"), "ImageJ-win64.exe")
	script_path = os.path.join(getProperty("fiji.dir"),  "plugins", "tribolium_processing", "tribolium_CLIJ_deconv_fusion.py")

	views_paths = ";".join(raw_transformed_paths)
	psfs_paths = ";".join(psf_paths)


	script_params = 'views_paths="%s",psfs_paths="%s",fused_output_path="%s",num_deconv_iter="%s"' % (views_paths, psfs_paths, fused_output_path, num_deconv_iter)
	print(script_params)
	command = [fiji_exe, '--ij2', '--headless', '--console', '--run', script_path, script_params.replace('\\', '/')]
	p = subprocess.Popen(command)
	logging_broadcast("Started CLIJ deconv. fusion process with PID: %s" % p.pid)

#						Pipeline helper functions

def save_tiff(image, path):
	if os.path.exists(path):
		os.remove(path)
	if COMPRESS_ON_SAVE == False:
		fs = FileSaver(image)
		fs.saveAsTiff(path)
	else:
		IJ.run(image, "Bio-Formats Exporter", "save=%s export compression=zlib" % path)

def get_directions_dirs_for_folder(input_dir, ndirections):
	direction_dirs = []
	for i in range(1, ndirections + 1):
			new_dir = os.path.join(input_dir, "DR" + str(i).zfill(4))
			direction_dirs.append(new_dir)
			if not os.path.exists(new_dir):
				mkpath(new_dir)

	return direction_dirs

def get_raw_images_dir(datasets_dir, dataset_id):
	dirs = [name for name in os.listdir(
		datasets_dir) if os.path.isdir(os.path.join(datasets_dir, name))]

	this_dataset_dir = fnmatch.filter(dirs, "DS%04d*" % dataset_id)
	if len(this_dataset_dir) > 1:
		error_msg = "	Error: there are multiple directories for the dataset with ID: %04d." % dataset_id
		logging.info(error_msg)
		print(error_msg)
		return None
	if len(this_dataset_dir) == 0:
		error_msg = "	Error: there are no directories for the dataset with ID: %04d." % dataset_id
		logging.info(error_msg)
		print(error_msg)
		return None
	raw_images_dir = os.path.join(
		datasets_dir, this_dataset_dir[0], RAW_IMAGES_DIR_NAME)
	if not os.path.isdir(raw_images_dir):
		error_msg = "	Error: there are no %s directoriy for the dataset with ID: %04d." % (
			RAW_IMAGES_DIR_NAME, dataset_id)
		logging.info(error_msg)
		print(error_msg)
		return None
	return raw_images_dir

def is_dataset_ID_input_valid(dataset_id):
	if not isinstance(dataset_id, int):
		return False
	if dataset_id in range(10000):
		return True
	return False

def is_specimen_input_valid(specimens_per_direction):
	if not isinstance(specimens_per_direction, tuple):
		return False
	for i in specimens_per_direction:
		if not isinstance(i, int):
			return False
		if i not in range(1000):
			return False
	return True

def is_timepoints_range_to_fuse_valid(tp_range):
	start = tp_range[0]
	end = tp_range[1]
	if not isinstance(start, int) or not isinstance(end, int):
		return False
	if end < start:
		return False
	if start < 1:
		return False
	return True

def metadata_file_check_for_errors(datasets_meta, datasets_dir):
	for dataset in datasets_meta["datasets"]:
		if "ID" not in dataset:
			print("Error while parsing .json file. Did not find a dataset ID for one of the datsets. Exiting.")
			exit(1)
		dataset_id = dataset["ID"]
		if "head_direction" not in dataset:
			print("Error while parsing .json file: no head_direction for the dataset with ID: \"%s\". Exiting." % dataset_id)
			print(EXAMPLE_JSON_METADATA_FILE)
			exit(1)
		specimen_directions_in_channels = []
		for chan in range(3):
			channel = "channel_%s" % chan
			if channel in dataset:
				if "specimens_for_directions_1234" not in dataset[channel]:
					print("Error while parsing .json file: no specimens_for_directions_1234 field in channel_%s for the dataset with ID: \"%s\". Exiting." % (channel, dataset_id))
					print(EXAMPLE_JSON_METADATA_FILE)
					exit(1)
				specimen_directions_in_channels.append(
					tuple(dataset[channel]["specimens_for_directions_1234"]))
		specimen_directions_in_channels = tuple(specimen_directions_in_channels)

		if not is_dataset_ID_input_valid(dataset_id):
			print("Error while parsing .json file: not a valid dataset ID: \"%s\". Dataset ID should be an integer from 0 to 9999. Exiting." % dataset_id)
			print(EXAMPLE_JSON_METADATA_FILE)
			exit(1)
		for specimens in specimen_directions_in_channels:
			if not is_specimen_input_valid(specimens):
				print("Error while parsing .json file: not a valid specimen list \"%s\" for the dataset with ID: \"%s\". Exiting." % (
					specimens, dataset_id))
				print(EXAMPLE_JSON_METADATA_FILE)
				exit(1)
		if not dataset["head_direction"] in ["right", "left"]:
			print("Error while parsing .json file: not a valid head_direction \"%s\" for the dataset with ID: \"%s\". Exiting." % (
				dataset["head_direction"], dataset_id))
			print(EXAMPLE_JSON_METADATA_FILE)
			exit(1)
		if "timepoints_range_to_fuse" in dataset.keys() and not is_timepoints_range_to_fuse_valid(dataset["timepoints_range_to_fuse"]):
			print("Error while parsing .json file: not valid timepoints_range_to_fuse \"%s\" for the dataset with ID: \"%s\". Exiting." % (
				dataset["head_direction"], dataset_id))
			print(EXAMPLE_JSON_METADATA_FILE)
			exit(1)
		if specimen_directions_in_channels == ():
			print("Error while parsing .json file: no channels found for the dataset with ID: \"%s\". Exiting." % dataset_id)
			print(EXAMPLE_JSON_METADATA_FILE)
			exit(1)
		raw_images_dir = get_raw_images_dir(datasets_dir, dataset_id)
		if not raw_images_dir:
			print("Could not find raw images dir: %s." % raw_images_dir)
			exit(1)
	logging.info("Checked metadata file. No errors found.")

def logging_broadcast(string):
	print(string)
	logging.info(string)


#					Pipeline main process functions

def move_files(raw_images_dir, specimen_directions_in_channels, dataset_id, dataset_name_prefix):
	"""Splits the embryo images by direction and puts in separate channel/direction folders. 
	Renames the files to conform to FredericFile format.

	Args:
		raw_images_dir (str): full path to mixed raw images files
		specimens_per_direction ( ((4,3,2,1),(5,6,8,7)) ): a tuple, each element has info for a channel with this index. 
		Info is a tuple with correspondance speciment number and directions, where directions are corresponding index in the tuple.
		dataset_id (int): ID number of the dataset
		dataset_name_prefix (str): prefix that the user specifies for the dataset

	Raises:
		Exception: metadata does not contain the specimen extracted from the file names
	"""
	specimens_info = {}
	ndirections = len(specimen_directions_in_channels[0])
	for channel, directions in enumerate(specimen_directions_in_channels, start=1):
		chan_dir = os.path.join(raw_images_dir, "CH%04d" % channel)
		if not os.path.exists(chan_dir):
			os.mkdir(chan_dir)
		direction_dirs = get_directions_dirs_for_folder(chan_dir, ndirections)
		for direct, specimen in enumerate(directions, start=1):
			specimens_info[specimen] = {"channel": channel, "direction": direct}

	for file_name in os.listdir(raw_images_dir):
		file_path = os.path.join(raw_images_dir, file_name)

		if os.path.isdir(file_path):
			continue
		if not file_name.endswith((".tif", ".TIF")):
			continue

		specimen = int(
			file_name[file_name.find("SPC0") + 4: file_name.find("SPC0") + 6])
		time_point = int(
			file_name[file_name.find("TL") + 2: file_name.find("TL") + 6]) + 1
		# channel = int(file_name[file_name.find("CHN") + 3: file_name.find("TL") + 5]) + 1 # Not used for old mDSLM images

		if specimen not in specimens_info:
			raise ValueError("In the metadata for the dataset: DS %04d there is no entry for the specimen: %i" % (
				dataset_id, specimen))

		image_channel = specimens_info[specimen]["channel"]
		embryo_direction = specimens_info[specimen]["direction"]
		new_file_name = "%sDS%04dTP%04dDR%04dCH%04dPL(ZS).tif" % (dataset_name_prefix,
															dataset_id,
															time_point,
															embryo_direction,
															image_channel)
		os.rename(file_path, os.path.join(
			direction_dirs[embryo_direction - 1], new_file_name))
		logging.info("New file \n%s\n Full path:\n%s\n Original name: \n%s\n Original path: \n%s\n" % (new_file_name,
																								 os.path.join(
																									 direction_dirs[embryo_direction - 1], new_file_name),
																								 file_name,
																								 file_path))

def b_branch_processing(dataset_metadata_obj):

	raw_images_dir = dataset_metadata_obj.raw_images_dir_path
	ndirections = dataset_metadata_obj.number_of_directions
	root_dataset_dir = dataset_metadata_obj.root_dir

	auto_cropbox_error = False

	# Setting channel=1 to Loop one time over all dimensions in channel 1 to determine the dataset_maximal_crop_box_width
	channel = 1
	chan_dir_name = "CH%04d" % channel
	raw_images_direction_dirs = get_directions_dirs_for_folder(
		os.path.join(raw_images_dir, chan_dir_name), ndirections)
	meta_dir = os.path.join(root_dataset_dir, METADATA_DIR_NAME)
	meta_d_dirs = get_directions_dirs_for_folder(os.path.join(meta_dir, chan_dir_name), ndirections)
	tstack_dataset_dirs = get_directions_dirs_for_folder(
				os.path.join(root_dataset_dir, TSTACKS_DIR_NAME, chan_dir_name), ndirections)
	raw_cropped_dirs = get_directions_dirs_for_folder(os.path.join(
		root_dataset_dir, RAW_CROPPED_DIR_NAME, chan_dir_name), ndirections)

	# Loop one time over all dimensions in channel 1 to determine the dataset_maximal_crop_box_width
	# We have to make sure that all directions in all channels have the same crop box dimensions.
	dataset_minimal_crop_box_dims = (0, 0)
	for raw_dir, tstack_dir, m_dir, direction in zip(raw_images_direction_dirs, tstack_dataset_dirs, meta_d_dirs, range(1, ndirections + 1)):
		logging.info(
			"\tCreating crop templates for all directions in CH0001 to ensure that all crop templates have the same dimensions.")
		logging.info("\tChannel: %s Direction: %s In the loop over directions. This iteration operating on the following directories:" % (
			channel, direction))
		logging.info("\n\t\t\t\t\t\t%s\n\t\t\t\t\t\t%s\n\t\t\t\t\t\t%s\n" %
						(raw_dir, tstack_dir, m_dir))

		tstack_backup_dir = os.path.join(tstack_dir, "uncropped_backup")
		if not os.path.exists(tstack_backup_dir):
			os.mkdir(tstack_backup_dir)
		logging.info("\tChannel: %s Direction: %s Creating a stack of max Z-projections from raw stack." %
						(channel, direction))
		mproj_stack_file_name = get_any_tiff_name_from_dir(raw_dir)
		mproj_stack_file_name.plane = "(ZM)"
		mproj_stack_file_name.time_point = "(TS)"
		mproj_stack_path = os.path.join(
			tstack_backup_dir, mproj_stack_file_name.get_name())
		if not os.path.exists(mproj_stack_path):
			max_proj_stack = make_max_Z_projection_stack_from_folder(raw_dir)
			fs = FileSaver(max_proj_stack)
			fs.saveAsTiff(mproj_stack_path)
		else:
			max_proj_stack = IJ.openImage(mproj_stack_path)
			logging.info(
				"\tChannel: %s Direction: %s Found existing max Z-projections. Using them." % (channel, direction))

		max_time_proj_file_name = mproj_stack_file_name
		max_time_proj_file_name.time_point = "(TM)"
		max_time_proj_full_path = os.path.join(
			m_dir, max_time_proj_file_name.get_name())
		if not os.path.exists(max_time_proj_full_path):
			max_time_proj = z_project_a_stack(max_proj_stack)
			fs = FileSaver(max_time_proj)
			fs.saveAsTiff(max_time_proj_full_path)
		else:
			logging.info(
				"\tChannel: %s Direction: %s Found existing max TIME-projection, using it. \n\t\t\t\t\t%s" % (channel, direction, max_time_proj_full_path))
			max_time_proj = IJ.openImage(max_time_proj_full_path)

		logging.info("\tChannel: %s Direction: %s Creating a crop template from a stack of max projections." % (
			channel, direction))
		if not auto_cropbox_error:
			try:
				crop_template, cropped_max_time_proj, dataset_minimal_crop_box_dims = create_crop_template(
					max_time_proj, m_dir, dataset_metadata_obj, dataset_minimal_crop_box_dims, use_dataset_box_dims=False)
			except Exception as e:
				logging_broadcast(
					"ERROR: Encountered an exception while trying to create a crop template. Skipping the dataset. Exception:\n%s\nContinuing to create max projections for the manual crop box." % e)
				auto_cropbox_error = True
				continue
			fs = FileSaver(cropped_max_time_proj)
			fs.saveAsTiff(os.path.join(
				m_dir, "cropped_max_time_proj.tif"))
	
	if auto_cropbox_error:
		return False

	# Main loop over the channels and directions for the dataset
	for channel in range(1, dataset_metadata_obj.number_of_channels + 1):
		chan_dir_name = "CH%04d" % channel
		raw_images_direction_dirs = get_directions_dirs_for_folder(
			os.path.join(raw_images_dir, chan_dir_name), ndirections)
		root_dataset_dir = os.path.split(raw_images_dir)[0]
		meta_dir = os.path.join(root_dataset_dir, METADATA_DIR_NAME)
		meta_d_dirs = get_directions_dirs_for_folder(os.path.join(meta_dir, chan_dir_name), ndirections)
		tstack_dataset_dirs = get_directions_dirs_for_folder(
			os.path.join(root_dataset_dir, TSTACKS_DIR_NAME, chan_dir_name), ndirections)
		raw_cropped_dirs = get_directions_dirs_for_folder(os.path.join(
			root_dataset_dir, RAW_CROPPED_DIR_NAME, chan_dir_name), ndirections)
		contrast_dirs = get_directions_dirs_for_folder(os.path.join(
			root_dataset_dir, CONTRAST_DIR_NAME, chan_dir_name), ndirections)

		for direction, raw_dir, tstack_dir, m_dir, raw_cropped_dir in zip(range(1, ndirections + 1), raw_images_direction_dirs, tstack_dataset_dirs, meta_d_dirs, raw_cropped_dirs):
			# if direction > 1:
			# 	print "Skipping all other directions for testing"
			# 	break
			logging.info("\tChannel: %s Direction: %s In the loop over directions. This iteration operating on the following directories:" % (
				channel, direction))
			logging.info("\n\t\t\t\t\t\t%s\n\t\t\t\t\t\t%s\n\t\t\t\t\t\t%s\n\t\t\t\t\t\t%s\n" % (
				raw_dir, tstack_dir, m_dir, raw_cropped_dir))
			tstack_backup_dir = os.path.join(tstack_dir, "uncropped_backup")

			if not os.path.exists(tstack_backup_dir):
				os.mkdir(tstack_backup_dir)
			logging.info(
				"\tChannel: %s Direction: %s Creating a stack of max Z-projections from raw stack." % (channel, direction))
			mproj_stack_file_name = get_any_tiff_name_from_dir(raw_dir)
			mproj_stack_file_name.plane = "(ZM)"
			mproj_stack_file_name.time_point = "(TS)"
			mproj_stack_path = os.path.join(
				tstack_backup_dir, mproj_stack_file_name.get_name())
			if not os.path.exists(mproj_stack_path):
				max_proj_stack = make_max_Z_projection_stack_from_folder(raw_dir)
				save_tiff(max_proj_stack, mproj_stack_path)
			else:
				max_proj_stack = IJ.openImage(mproj_stack_path)
				logging.info(
					"\tChannel: %s Direction: %s Found existing max Z-projections. Using them." % (channel, direction))

			max_time_proj_file_name = mproj_stack_file_name
			max_time_proj_file_name.time_point = "(TM)"
			max_time_proj_full_path = os.path.join(
				m_dir, max_time_proj_file_name.get_name())
			if not os.path.exists(max_time_proj_full_path):
				max_time_proj = z_project_a_stack(max_proj_stack)
				fs = FileSaver(max_time_proj)
				fs.saveAsTiff(max_time_proj_full_path)
			else:
				logging.info(
					"\tChannel: %s Direction: %s Found existing max TIME-projection, using it. \n\t\t\t\t\t%s" % (channel, direction, max_time_proj_full_path))
				max_time_proj = IJ.openImage(max_time_proj_full_path)

			logging.info("\tChannel: %s Direction: %s Creating a crop template from a stack of max projections." % (
				channel, direction))
			try:
				crop_template, cropped_max_time_proj, dataset_minimal_crop_box_dims = create_crop_template(
					max_time_proj, m_dir, dataset_metadata_obj, dataset_minimal_crop_box_dims, use_dataset_box_dims=True)
			except Exception as e:
				traceback.print_exc()
				logging.info(
					"ERROR: Encountered an exception while trying to create a crop template. Skipping the dataset. Exception:\n%s" % e)
				print("ERROR: Encountered an exception while trying to create a crop template. Skipping the dataset. Exception:\n%s" % e)
				return False
			fs = FileSaver(cropped_max_time_proj)
			fs.saveAsTiff(os.path.join(
				m_dir, "cropped_max_time_proj.tif"))

			# Cropping max projections

			cropped_tstack_file_name = get_any_tiff_name_from_dir(
				raw_dir)
			cropped_tstack_file_name.time_point = "(TS)"
			cropped_tstack_file_name.plane = "(ZM)"
			if not USE_CROPPED_CACHE or not os.path.exists(os.path.join(tstack_dir, cropped_tstack_file_name.get_name())):
				logging.info("\tChannel: %s Direction: %s Cropping a stack of max projections." % (
					channel, direction))
				cropped_tstack = crop_stack_by_template(
					max_proj_stack, crop_template, dataset_metadata_obj)
				save_tiff(cropped_tstack, os.path.join(
					tstack_dir, cropped_tstack_file_name.get_name()))
			else:
				logging.info("\tChannel: %s Direction: %s Found existing cropped max projected stacks. Using them." % (
											channel, direction))
				cropped_tstack = IJ.openImage(os.path.join(
										tstack_dir, cropped_tstack_file_name.get_name()))

			# Cropping raw stacks

			ntime_points = len(get_tiffs_in_directory(raw_dir))
			stacks_done_list = get_tiffs_in_directory(raw_cropped_dir)
			all_stacks_done = ntime_points == len(stacks_done_list)
			if USE_CROPPED_CACHE and all_stacks_done:
				logging.info("\tChannel: %s Direction: %s Found existing cropped raw stacks, using them." % (
										channel, direction))
			else:
				logging.info(
					"\tChannel: %s Direction: %s Cropping raw image stacks." % (channel, direction))
				planes_kept = (0, 0)
				raw_stack_files = get_tiffs_in_directory(raw_dir)
				if len(stacks_done_list) > 1:
					logging.info("\tChannel: %s Direction: %s Found some existing cropped raw stacks, using them." % (
												channel, direction))
					files_to_crop = raw_stack_files[len(stacks_done_list) - 1 : ]
					current_active_stack = len(stacks_done_list)
				else:
					files_to_crop = raw_stack_files
					current_active_stack = 1
				for direction, raw_stack_file_name in enumerate(files_to_crop):
					raw_stack = IJ.openImage(raw_stack_file_name)
					IJ.run(raw_stack, "Properties...",
											"frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=4.0000")
					raw_stack_cropped = crop_stack_by_template(
						raw_stack, crop_template, dataset_metadata_obj)
					if direction == 0:
						logging.info("\tChannel: %s Direction: %s Finding which planes to keep in raw image stacks." % (
							channel, direction))
						if dataset_metadata_obj.planes_to_keep_per_direction_b_branch is not None:
							manual_planes_to_keep = (dataset_metadata_obj.planes_to_keep_per_direction_b_branch[direction - 1]["start"],
							 dataset_metadata_obj.planes_to_keep_per_direction_b_branch[direction - 1]["end"])
						else:
							manual_planes_to_keep = None
						stack_for_finding_planes = IJ.openImage(raw_stack_files[0])
						IJ.run(stack_for_finding_planes, "Properties...",
											"frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=4.0000")
						stack_for_finding_planes = crop_stack_by_template(
														stack_for_finding_planes, crop_template, dataset_metadata_obj)
						try:
							planes_kept = find_planes_to_keep(stack_for_finding_planes, m_dir, manual_planes_to_keep)
						except Exception as e:
							logging.info("\tChannel: %s Direction: %s Encountered an exception while trying to find which planes to keep. Skipping the dataset. Exception: \n%s" % (
								channel, direction, e))
							print("\tChannel: %s Direction: %s Encountered an exception while trying to find which planes to keep. Skipping the dataset. Exception: \n%s" % (
								channel, direction, e))
							print("Encountered an exception while trying to find which planes to keep. Skipping the dataset B-branch.")
							return False
						logging.info("\tChannel: %s Direction: %s Keeping planes: %s." %
									(channel, direction, planes_kept))
					logging.info(
						"\tChannel: %s Direction: %s Cropping Raw stack for timepoint: %s/%s" % (channel, direction, current_active_stack, ntime_points))
					raw_stack_cropped = reset_img_properties(raw_stack_cropped, voxel_depth=4)
					raw_stack_cropped = subset_planes(raw_stack_cropped, planes_kept)
					save_tiff(raw_stack_cropped, os.path.join(
						raw_cropped_dir, os.path.split(raw_stack_file_name)[1]))
					current_active_stack += 1


		# Save unadjusted montages

		logging.info("\tChannel: %s Creating montages from max projections." % (
						channel))
		montage_stack = ImagePlus()
		montage_stack_name = None
		for direction, tstack_dir in enumerate(tstack_dataset_dirs):
			stack_path = os.path.join(
				tstack_dir, get_any_tiff_name_from_dir(tstack_dir).get_name())
			imp_stack = IJ.openImage(stack_path)
			if direction == 0:
				montage_stack = imp_stack.getStack()
				montage_stack_name = get_any_tiff_name_from_dir(tstack_dir)
			if direction > 0:
				montage_stack = StackCombiner.combineHorizontally(StackCombiner(), montage_stack, imp_stack.getStack())
		montage_stack_name.direction = "(AX)"
		montage_stack_name.plane = "(ZM)"
		montage_stack = ImagePlus("montage", montage_stack)
		montage_dir = os.path.join(root_dataset_dir, MONTAGE_DIR_NAME)
		if not os.path.exists(montage_dir):
			mkpath(montage_dir)
		save_tiff(montage_stack, os.path.join(montage_dir, montage_stack_name.get_name()))


		# Creating bleach corrected and adjusted montages
			
		logging.info("\tChannel: %s Creating histogram adjusted montages from max projections." % (
									channel))
		adj_montage_stack = montage_stack
		adj_montage_stack_name = montage_stack_name
		adj_montage_stack_name.direction = "(AX)"
		adj_montage_stack_name.plane = "(ZA)"
		if DO_HISTOGRAM_MATCHING:
			match_histograms_stack(adj_montage_stack)
			adj_montage_stack_name.plane = "(ZH)"
		adj_montage_stack, histogram_thresholds = threshold_histogram_stack(adj_montage_stack, PERCENT_OVEREXPOSED_PIXELS)
		logging.info("Adjusted stack with thresholds: %s, %s." % (
			histogram_thresholds[0], histogram_thresholds[1]))
		histogram_thresholds = {"lower_threshold": histogram_thresholds[0], "upper_threshold": histogram_thresholds[1]}
		with open(os.path.join(meta_dir, chan_dir_name, "histogram_thresholds_used_for_contrast_adjustment.JSON"), "w") as hist_json:
			json.dump(histogram_thresholds, hist_json, indent=4)
		
		save_tiff(adj_montage_stack, os.path.join(montage_dir, adj_montage_stack_name.get_name()))


		# Split the adjusted montage into separate stacks for directions
		logging.info("\tChannel: %s Saving histogram adjusted stacks of max projections." % (
									channel))
		adjusted_stacks = split_montage_stack_to_list_of_stacks(adj_montage_stack, 1, ndirections)

		for direction, tstack_dir, contr_dir, adj_max_proj_stack in zip(range(1, ndirections + 1), tstack_dataset_dirs, contrast_dirs, adjusted_stacks):
			logging.info("\tChannel: %s Direction: %s Saving histogram adjusted stack." % (
														channel, direction))
			adj_stack_name = get_any_tiff_name_from_dir(tstack_dir)
			adj_stack_name.plane = "(ZA)"
			if DO_HISTOGRAM_MATCHING:
				adj_stack_name.plane = "(ZH)"
		
			save_tiff(adj_max_proj_stack, os.path.join(contr_dir, adj_stack_name.get_name()))
		
		# Save adjusted vertical montages

		logging.info("\tChannel: %s Creating adjusted vertical montages from max projections." % (
						channel))
		vert_montage_stack = ImagePlus()
		vert_montage_stack_name = None
		for direction, tstack_dir in enumerate(tstack_dataset_dirs):
			stack_path = os.path.join(
				tstack_dir, get_any_tiff_name_from_dir(tstack_dir).get_name())
			imp_stack = IJ.openImage(stack_path)
			if direction == 0:
				vert_montage_stack = imp_stack.getStack()
				vert_montage_stack_name = get_any_tiff_name_from_dir(tstack_dir)
			if direction > 0:
				vert_montage_stack = StackCombiner.combineVertically(StackCombiner(), vert_montage_stack, imp_stack.getStack())
		vert_montage_stack_name.direction = "(AY)"
		vert_montage_stack_name.plane = "(ZA)"
		vert_montage_stack = ImagePlus("montage", vert_montage_stack)
		if DO_HISTOGRAM_MATCHING:
			match_histograms_stack(vert_montage_stack)
			vert_montage_stack_name.plane = "(ZH)"
		vert_montage_stack, histogram_thresholds = threshold_histogram_stack(vert_montage_stack, PERCENT_OVEREXPOSED_PIXELS)
		logging.info("Adjusted stack with thresholds: %s, %s." % (
			histogram_thresholds[0], histogram_thresholds[1]))
		montage_dir = os.path.join(root_dataset_dir, MONTAGE_DIR_NAME)
		if not os.path.exists(montage_dir):
			mkpath(montage_dir)
		save_tiff(vert_montage_stack, os.path.join(montage_dir, vert_montage_stack_name.get_name()))



	return True

def fusion_branch_processing(dataset_metadata_obj):
	raw_img_ch_1_dir = os.path.join(dataset_metadata_obj.raw_images_dir_path, "CH0001")
	raw_direction_dirs = get_directions_dirs_for_folder(raw_img_ch_1_dir, dataset_metadata_obj.number_of_directions)
	example_raw_filename = get_any_tiff_name_from_dir(raw_direction_dirs[0])
	reference_timepoint = dataset_metadata_obj.tp_fuse_start


	downsampled_dir = os.path.join(dataset_metadata_obj.root_dir, "downsampled_%s_timepoint" % reference_timepoint)
	if not os.path.exists(downsampled_dir):
		mkpath(downsampled_dir)

	cropped_projections_dir = os.path.join(dataset_metadata_obj.datasets_dir, "cropped_fused_projections")
	if not os.path.exists(cropped_projections_dir):
		mkpath(cropped_projections_dir)

	view_image_filename = example_raw_filename
	view_image_filename.time_point = "%04d" % reference_timepoint

	# Downsample reference timepoint to a new folder
	logging_broadcast("Downsampling images for fusion+segmentation")
	for direction, direction_dir in zip(range(1, dataset_metadata_obj.number_of_directions + 1), raw_direction_dirs):
		view_image_filename.direction = "%04d" % (direction + 1)
		next_downsampled_image_path = os.path.join(downsampled_dir, view_image_filename.get_name())

		view_image_filename.direction = "%04d" % direction
		downsampled_image_path = os.path.join(downsampled_dir, view_image_filename.get_name())

		if os.path.exists(downsampled_image_path) and os.path.exists(next_downsampled_image_path): 
			# We cannot guarantee that last downsampled file was finished if the pipeline crashed previously, so always redo it
			logging_broadcast("Found existing file %s, skipping downsampling." % downsampled_image_path)
			continue

		full_image_path = os.path.join(direction_dir, view_image_filename.get_name())
		full_image = IJ.openImage(full_image_path)
		downsampled_image = downsample_image_stack(full_image, full_image.getWidth() / IMAGE_DOWNSAMPLING_FACTOR_XY)
		fs = FileSaver(downsampled_image)
		fs.saveAsTiff(downsampled_image_path)

	# Fuse downsampled reference timepoint
	view_dims = get_image_dimensions_from_file(os.path.join(downsampled_dir, view_image_filename.get_name()))
	fuse_file_pattern = view_image_filename
	fuse_file_pattern.direction = "{aaaa}"
	logging_broadcast("DS%04d Fusing downsampled %s timepoint as reference" % (dataset_metadata_obj.id, reference_timepoint))
	if not DEBUG_USE_FUSED_PREVIEW_CACHE:
		create_and_fuse_preview_dataset(downsampled_dir, fuse_file_pattern.get_name(), (1, dataset_metadata_obj.number_of_directions), view_dims)

	fused_file_path = os.path.join(downsampled_dir, "fused_tp_0_ch_0.tif")

	fused_stack = IJ.openImage(fused_file_path)
	z_projection = tribolium_image_utils.project_image(fused_stack, "Z", "Max")
	y_projection = tribolium_image_utils.project_image(fused_stack, "Y", "Max")

	logging_broadcast("Segmenting downsampled and fusing cropped embryo")
	zy_cropped_projections, segmentation_results = segment_embryo_and_fuse_again_cropping_around_embryo(
														raw_dataset_xml_path=os.path.join(downsampled_dir, "dataset.xml"),
														fused_xml_path=os.path.join(downsampled_dir, "fused.xml"),
														z_projection=z_projection,
														y_projection=y_projection)
	
	# Save ZY cropped projection montages
	if zy_cropped_projections == False:
		logging_broadcast("segment_embryo_and_fuse_again_cropping_around_embryo Failed. Exiting Fusion-branch.")
		return False
	IJ.setBackgroundColor(255, 255, 255)
	IJ.run(zy_cropped_projections, "Canvas Size...", "width=%s height=%s position=Center" % (CANVAS_SIZE_FOR_EMBRYO_PREVIEW, CANVAS_SIZE_FOR_EMBRYO_PREVIEW))
	zy_projections_filename = view_image_filename
	zy_projections_filename.direction = "z+yM"
	zy_projections_external_path = os.path.join(cropped_projections_dir, zy_projections_filename.get_name_without_extension() + "_montage.tiff")
	zy_projections_meta_path = os.path.join(dataset_metadata_obj.root_dir, METADATA_DIR_NAME, zy_projections_filename.get_name_without_extension() + "_montage.tiff")
	save_tiff(zy_cropped_projections, zy_projections_external_path)
	save_tiff(zy_cropped_projections, zy_projections_meta_path)

	if RUN_FUSION_UP_TO_DOWNSAMPLED_PREVIEW:
		logging_broadcast("Fusion-branch was specified to run only until preview. Continuing.")
		return True

	###### Fusion and deconvolution full dataset

	upscaled_translation_embryo_to_center = upscale_translation_matrix(segmentation_results.transformation_embryo_to_center, IMAGE_DOWNSAMPLING_FACTOR_XY)
	upscaled_bounding_box = upscale_bounding_box(segmentation_results.embryo_crop_box, IMAGE_DOWNSAMPLING_FACTOR_XY, (IMAGE_PIXEL_DISTANCE_X, IMAGE_PIXEL_DISTANCE_Y, IMAGE_PIXEL_DISTANCE_Z))
	upscaled_vertical_box = rotate_bounding_box_vertically(upscaled_bounding_box)
	
	logging.info("Upscaled matrix for embryo translation to center: %s" % upscaled_translation_embryo_to_center)
	logging.info("Upscaled bounding box for embryo cropping: %s" % upscaled_bounding_box)

	for ch in range(1, dataset_metadata_obj.number_of_channels + 1):
		raw_direction_dirs = get_directions_dirs_for_folder(os.path.join(dataset_metadata_obj.raw_images_dir_path, "CH%04d" % ch), 
															dataset_metadata_obj.number_of_directions)
		fusion_setup_folder = os.path.join(dataset_metadata_obj.root_dir, RAW_IMAGES_DIR_NAME, "CH%04d" % ch, "bigstitcher_dataset")
		if not os.path.exists(fusion_setup_folder):
			mkpath(fusion_setup_folder)
		if not os.path.exists(dataset_metadata_obj.raw_aligned_stacks_dir):
			mkpath(dataset_metadata_obj.raw_aligned_stacks_dir)


		pre_fusion_dataset_basename = view_image_filename
		pre_fusion_dataset_basename.direction = "all_"
		pre_fusion_dataset_basename.time_point = "all_"
		pre_fusion_dataset_basename.channel = "%04d" % ch
		fusion_dataset_basename = deepcopy(pre_fusion_dataset_basename)
		fusion_dataset_basename.direction = "fuse"
		fusion_dataset_basename = os.path.splitext(fusion_dataset_basename.get_name())[0]
		full_dataset_file_pattern = deepcopy(pre_fusion_dataset_basename)
		pre_fusion_dataset_basename = os.path.splitext(pre_fusion_dataset_basename.get_name())[0]
		full_dataset_file_pattern.direction = "{aaaa}"
		full_dataset_file_pattern.time_point = "{tttt}"
		full_dataset_file_pattern = os.path.join("..", "DR{aaaa}", full_dataset_file_pattern.get_name())

		selected_timepoints = get_timepoint_list_from_folder(
											raw_direction_dirs[0], 
											dataset_metadata_obj.tp_fuse_start, 
											dataset_metadata_obj.tp_fuse_end)
		
		dataset_xml_path = os.path.join(fusion_setup_folder, pre_fusion_dataset_basename + ".xml")
		is_bounding_box_present = get_bounding_box_coords_from_xml(dataset_xml_path, "embryo_cropped") is not None
		
		if USE_FUSION_DATASET_CACHE and is_bounding_box_present:
			logging_broadcast("Found previously created dataset with 'embryo_cropped' crop box. Skipping dataset creation and registration.")
		else:
			logging_broadcast("Creating dataset for fusion with selected timepoints: %s-%s" % (selected_timepoints[0], selected_timepoints[-1]))
			create_and_register_full_dataset(fusion_setup_folder, 
											pre_fusion_dataset_basename, 
											full_dataset_file_pattern, 
											selected_timepoints,
											(1, dataset_metadata_obj.number_of_directions),
											(IMAGE_PIXEL_DISTANCE_X, IMAGE_PIXEL_DISTANCE_Y, IMAGE_PIXEL_DISTANCE_Z),
											reference_timepoint=reference_timepoint)
			apply_transformation_bigstitcher_dataset(dataset_xml_path, upscaled_translation_embryo_to_center)
			rotate_bigstitcher_dataset(dataset_xml_path, "z", segmentation_results.rotation_angle_z)
			rotate_bigstitcher_dataset(dataset_xml_path, "y", segmentation_results.rotation_angle_y)
			if dataset_metadata_obj.embryo_head_direction == "right":
				rotate_bigstitcher_dataset(dataset_xml_path, "z", -90)
			if dataset_metadata_obj.embryo_head_direction == "left":
				rotate_bigstitcher_dataset(dataset_xml_path, "z", 90)
			define_bounding_box_for_fusion(dataset_xml_path, upscaled_vertical_box, "embryo_cropped")


		print("Starting Fusion or Deconvolution")
		fusion_output_dir = os.path.join(dataset_metadata_obj.root_dir, FUSED_DIR_NAME, "CH%04d" % ch)
		if not os.path.exists(fusion_output_dir):
			mkpath(fusion_output_dir)
		if RUN_FUSION_BRANCH_UP_TO_START_OF_FUSION:
			logging_broadcast("Specified to run only up to the start of deconvolution. Done fusion branch.")
			return True
		if ONLY_FUSE_FULL:
			fuse_dataset_to_tiff(dataset_xml_path, "embryo_cropped", os.path.join(fusion_output_dir, fusion_dataset_basename + ".xml"))
		if FUSE_CONTENT_BASED:
			fuse_dataset_to_tiff(dataset_xml_path, "embryo_cropped", os.path.join(fusion_output_dir, fusion_dataset_basename + ".xml"), use_entropy_weighted_fusion=True)
		if FUSE_AND_DECONVOLVE_FULL:
			extract_psf(dataset_xml_path, reference_timepoint)
			assign_psf(dataset_xml_path, reference_timepoint)
			deconvolve_dataset_to_tiff(dataset_xml_path, 
										"embryo_cropped",
										os.path.join(fusion_output_dir, fusion_dataset_basename + ".xml"),
										NUM_DECONV_ITERATIONS)
		if FUSE_CLIJ_GPU:
			if USE_FUSION_DATASET_CACHE and os.path.exists(os.path.join(fusion_setup_folder, "psf")):
				logging_broadcast("Found 'psf' folder inside dataset setup folder. Skipping PSF extraction and asignment.")
			else:	
				logging_broadcast("Extracting PSF from reference timepoint %s and assigning it to all." % reference_timepoint)
				extract_psf(dataset_xml_path, reference_timepoint)
				assign_psf(dataset_xml_path, reference_timepoint)

			# temp_dir_fusion = os.path.join(CACHING_DIR, uuid.uuid4().hex) # CHANGE BACK!!!! #####################################
			# mkpath(temp_dir_fusion)
			temp_dir_fusion = os.path.join(CACHING_DIR, "fusion_cache_dir_DS00" + str(dataset_metadata_obj.id) + "_" + os.path.basename(dataset_metadata_obj.datasets_dir)) #  uuid.uuid4().hex # CHANGE BACK!!!! #####################################
			if not os.path.exists(temp_dir_fusion):
				mkpath(temp_dir_fusion)	

			logging_broadcast("Starting deconvolution->content-based fusion on the GPU ")
			initial_mem_usage = get_free_memory_in_GB()
			used_mem_last_iter = initial_mem_usage
			last_fused_path = ""
			last_raw_transformed_paths = []
			for tp_index, tp in enumerate(selected_timepoints):
				fused_save_path =  os.path.join(fusion_output_dir, "deconv_weighted_fused_tp_%s.tiff" % tp)
				logging_broadcast("Fusing timepoint %s" % tp)
				if os.path.exists(fused_save_path):
					logging_broadcast("Found previously fused timepoint, skipping.")
					last_fused_path = "skipped"
					continue

				transformed_psf_dir = os.path.join(fusion_output_dir, "transformed_psf")
				if not os.path.exists(transformed_psf_dir):
					mkpath(transformed_psf_dir)
				logging_broadcast("Saving averaged transformed PSFs")
				psf_paths = save_transformed_psfs(dataset_xml_path, transformed_psf_dir, tp, dataset_metadata_obj.number_of_directions)


				logging_broadcast("Saving raw transformed stacks")
				raw_transformed_paths = save_raw_transformed_stacks(dataset_xml_path, temp_dir_fusion, tp, dataset_metadata_obj.number_of_directions)
				if tp_index > 0:
					logging_broadcast("Checking wheather last timepoint has fused successfully. Expected path: %s  Is file present: %s " % (last_fused_path, os.path.exists(last_fused_path)))
					if last_fused_path != "skipped" and not os.path.exists(last_fused_path):
						logging_broadcast("ERROR: Prevous timepoint fusion was not generated for some reason. Exiting pipeline!")
						exit(1)
					if os.path.exists(last_fused_path):
						logging_broadcast("Started moving raw aligned stack for previous fused timepoint.")
						start_moving_files(last_raw_transformed_paths, dataset_metadata_obj.raw_aligned_stacks_dir)


				logging_broadcast("Starting deconvolution->fusion of the timepoint with params raw_transformed_paths: %s psf_paths: %s fused_save_path: %s NUM_DECONV_ITERATIONS: %s" % (raw_transformed_paths, psf_paths, fused_save_path, NUM_DECONV_ITERATIONS))
				# Spawning a separate process for fusion and deconvolution because this function leaks memory, and this is the only way to make sure it is released.
				start_deconv_fuse_timepoint_process(raw_transformed_paths, psf_paths, fused_save_path, NUM_DECONV_ITERATIONS)

				# shutil.rmtree(temp_dir_fusion) ## Change this! ##########################################
				used_mem = get_free_memory_in_GB()
				logging_broadcast("Used memory growth per iteration: %s Gb" % (used_mem - used_mem_last_iter))
				logging_broadcast("Used memory growth total: %s Gb" % (used_mem - initial_mem_usage))
				logging_broadcast("Used memory total: %s Gb" % used_mem)
				used_mem_last_iter = used_mem
				last_fused_path = fused_save_path
				last_raw_transformed_paths = raw_transformed_paths
			for sec in range(400):
				if os.path.exists(last_fused_path):
					logging_broadcast("Started moving raw aligned stack for the last timepoint.")
					start_moving_files(last_raw_transformed_paths, dataset_metadata_obj.raw_aligned_stacks_dir)
					return True
				time.sleep(1)
			logging_broadcast("ERROR: reached timeout on waiting for the last fused timepoint.")
			return False
				

	return True

def process_datasets(datasets_dir, metadata_file, dataset_name_prefix):
	"""Main function for processing the datasets.

	Args:
		datasets_dir (java.io.File): directory with subdirectories for each dataset	(input from the user)
		metadata_file (java.io.File): JSON file with metadata for each dataset (input from the user)
		dataset_name_prefix (str): Prefix that will be added to all image files for all datasets (input from the user)
	"""
	# Converting a Java File object to a string.
	if isinstance(metadata_file, File):
		metadata_file = metadata_file.getAbsolutePath()
	if isinstance(datasets_dir, File):
		datasets_dir = datasets_dir.getAbsolutePath()

	with open(metadata_file) as f:
		try:
			datasets_meta = json.load(f)
		except ValueError as err:
			print("Could not load the JSON metadata file.")
			print("Error generated by the JSON parser: \n%s" % err)
			print(EXAMPLE_JSON_METADATA_FILE)
			exit(1)

	now = datetime.now()
	dt_string = now.strftime("%Y-%b-%d-%H%M%S")
	logging.basicConfig(filename=os.path.join(datasets_dir, "%s-TribFuse_pipeline.log" % dt_string),
					 filemode='w',
					 format='%(asctime)s-%(levelname)s - %(message)s',
					 datefmt='%d-%b-%y %H:%M:%S',
					 level=logging.INFO)

	metadata_file_check_for_errors(datasets_meta=datasets_meta, datasets_dir=datasets_dir)
	if not os.path.exists(CACHING_DIR) and CACHING_DIR != None:
		logging_broadcast("Could not find cache directory specified. Exiting. Directory provided: %s" % CACHING_DIR)
		exit(1)



	for dataset in datasets_meta["datasets"]:
		skip_the_dataset = False

		dataset_metadata_obj = get_DatasetMetadata_obj_from_metadata_dict(dataset, datasets_dir)
		raw_images_dir = dataset_metadata_obj.raw_images_dir_path
		root_dataset_dir = dataset_metadata_obj.root_dir

		specimen_directions_in_channels = dataset_metadata_obj.specimen_directions_in_channels
		dataset_id = dataset_metadata_obj.id

		skip_b_branch = False
		skip_fusion_branch = False

		
		logging_broadcast("\n%s\nStarted processing dataset: DS%04d \n%s\n" %
					 ("#" * 100, dataset_id, "#" * 100))
		if os.path.exists(os.path.join(root_dataset_dir, B_BRANCH_FINISHED_FILE_NAME)):
			logging_broadcast(
				"Found %s file. Dataset DS%04d B-branch already processed, skipping." % (B_BRANCH_FINISHED_FILE_NAME, dataset_id))
			skip_b_branch = True
		if os.path.exists(os.path.join(root_dataset_dir, FUSION_BRANCH_FINISHED_FILE_NAME)):
			logging_broadcast(
				"Found %s file. Dataset DS%04d FUSION-branch already processed, skipping." % (FUSION_BRANCH_FINISHED_FILE_NAME, dataset_id))
			skip_fusion_branch = True

		if os.path.exists(os.path.join(root_dataset_dir, B_BRANCH_ERROR_FILE_NAME)):
			logging_broadcast(
				"Found %s file. Dataset DS%04d B-branch errored while previous processing, skipping." % (B_BRANCH_ERROR_FILE_NAME, dataset_id))
			skip_b_branch = True
		if os.path.exists(os.path.join(root_dataset_dir, FUSION_BRANCH_ERROR_FILE_NAME)):
			logging_broadcast(
				"Found %s file. Dataset DS%04d FUSION-branch errored while previous processing, skipping." % (FUSION_BRANCH_ERROR_FILE_NAME, dataset_id))
			skip_fusion_branch = True


		# Here there maybe a race condition for moving files if two pipeline instances are spawned in parallel. But it should not be a problem.
		logging_broadcast("\tArranging raw image files")
		try:
			move_files(
				raw_images_dir, specimen_directions_in_channels, dataset_id, dataset_name_prefix)
		except ValueError as e:
			logging_broadcast(
				"Error while moving files for the dataset:\"%s\", skipping the dataset. Error:\n %s" % (dataset_id, e))
			continue
		logging_broadcast("\tFiles arranged, starting processing.")

		
		if os.path.exists(os.path.join(root_dataset_dir, B_BRANCH_LOCK_FILE_NAME)):
			logging_broadcast(
				"Found %s file. Perhaps the dataset DS%04d B-branch is currently being processed by other Fiji instance, skipping." % (B_BRANCH_LOCK_FILE_NAME, dataset_id))
			skip_b_branch = True
		if DO_B_BRANCH and not skip_b_branch:
			open(os.path.join(root_dataset_dir, B_BRANCH_LOCK_FILE_NAME), 'a').close()
			logging_broadcast("Starting B-branch processing. DS%04d" % dataset_metadata_obj.id)
			b_branch_successfull = b_branch_processing(dataset_metadata_obj=dataset_metadata_obj)
			os.remove(os.path.join(root_dataset_dir, B_BRANCH_LOCK_FILE_NAME))
			if not b_branch_successfull:
				logging_broadcast("ERROR: B-branch failed. Had to skip the dataset DS%04d." % dataset_metadata_obj.id)
				open(os.path.join(root_dataset_dir, B_BRANCH_ERROR_FILE_NAME), 'a').close()
			else:
				open(os.path.join(root_dataset_dir, B_BRANCH_FINISHED_FILE_NAME), 'a').close()
		if not DO_B_BRANCH: 
			logging_broadcast("B-branch processing set to FALSE. Skipping B-branch. DS%04d." % dataset_metadata_obj.id)



		if os.path.exists(os.path.join(root_dataset_dir, FUSION_BRANCH_LOCK_FILE_NAME)):
			logging_broadcast(
				"Found %s file. Perhaps the dataset DS%04d FUSION-branch is currently being processed by other Fiji instance, skipping." % (FUSION_BRANCH_LOCK_FILE_NAME, dataset_id))
			skip_fusion_branch = True
		if DO_FUSION_BRANCH and not skip_fusion_branch:
			open(os.path.join(root_dataset_dir, FUSION_BRANCH_LOCK_FILE_NAME), 'a').close()
			logging_broadcast("Starting Fusion-branch processing. DS%04d" % dataset_metadata_obj.id)
			fusion_branch_successfull = fusion_branch_processing(dataset_metadata_obj=dataset_metadata_obj)
			os.remove(os.path.join(root_dataset_dir, FUSION_BRANCH_LOCK_FILE_NAME))
			if not fusion_branch_successfull:
				logging.info("ERROR: Fusion-branch failed. Had to skip the dataset DS%04d." % dataset_metadata_obj.id)
				open(os.path.join(root_dataset_dir, FUSION_BRANCH_ERROR_FILE_NAME), 'a').close()
			else:
				open(os.path.join(root_dataset_dir, FUSION_BRANCH_FINISHED_FILE_NAME), 'a').close()
		if not DO_FUSION_BRANCH: 
			logging_broadcast("Fusion-branch processing set to FALSE. Skipping Fusion-branch. DS%04d." % dataset_metadata_obj.id)

		logging_broadcast("Finished processing dataset DS%04d." % dataset_id)


	logging_broadcast("Finished processing all datasets.")


if __name__ in ['__builtin__', '__main__']:
	if DEBUG_USE_FUSED_PREVIEW_CACHE:
		print("DEBUGGING IS ACTIVE!!!")
	process_datasets(datasets_dir, metadata_file, dataset_name_prefix)

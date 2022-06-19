import math
import sys
import os

from java.lang.System import getProperty

from net.imagej.axis import Axes
from net.imagej.ops import Ops
from net.imagej import Dataset

from ij import IJ, ImagePlus, ImageStack
from fiji.threshold import Auto_Threshold
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from ij.plugin import ZProjector
from ij.plugin import StackCombiner
from ij.plugin import StackMaker
from ij.gui import RotatedRectRoi
from emblcmci import BleachCorrection_MH
from ij.plugin.filter import Analyzer
from ij.plugin import RoiEnlarger

opsServiceImageUtilsLocal = None
dsServiceImageUtilsLocal = None
convertServiceImageUtilsLocal = None

sys.path.append(os.path.join(getProperty("fiji.dir"), "plugins", "tribolium_processing")) 
from tribolium_math_utils import *
from tribolium_file_utils import *

###### Image managing

def z_project_a_stack(stack):
	"""Project a stack of images along the Z axis 

	Args:
		stack ImagePlus: stack of images to project

	Returns:
		ImagePlus: max-projection
	"""
	zp = ZProjector(stack)
	zp.setMethod(ZProjector.MAX_METHOD)
	zp.doProjection()
	zpimp = zp.getProjection()
	return zpimp

def is_polygon_roi_overlapping_image_edges(image, roi):
	"""Check if any of the roi coordinates are outside of the image.

	Args:
		image (ImagePlus): image to check over
		roi (PolygonRoi): roi to check

	Returns:
		bool: is roi overlapping the edges?
	"""
	is_overlapping = False
	if not all(x < image.getWidth() for x in roi.getPolygon().xpoints):
		is_overlapping = True
	if not all(x > 0 for x in roi.getPolygon().xpoints):
		is_overlapping = True
	if not all(y < image.getHeight() for y in roi.getPolygon().ypoints):
		is_overlapping = True
	if not all(y > 0 for y in roi.getPolygon().ypoints):
		is_overlapping = True
	return is_overlapping

def crop_stack_by_template(stack, crop_template, dataset_metadata_obj):
	"""Crop stack by provided PolygonRoi which should be a rotated rectangle around the embryo.
	It also flips the ebryo accroding to metadata about embryo facing direction, so that illumination is from the left.

	Args:
		stack (ImagePlus): a stack of image planes
		crop_template (PolygonRoi): rotated bounding rectangle around the embryo
		dataset_metadata_obj (DatasetMetadata): metadata about the dataset_metadata_obj to extract information about embryo orientation

	Returns:
		ImagePlus: cropped stack
	"""

	stack.setRoi(crop_template)
	cropped_stack = stack.crop("stack")
	stack = None
	# cropped_stack.show()
	IJ.run(cropped_stack, "Select All", "")

	IJ.run(cropped_stack, "Rotate... ",
			"angle=%s grid=1 interpolation=Bilinear stack" % round(get_rectangular_polygon_roi_angle(crop_template), ndigits=1))
	final_center_x = cropped_stack.getWidth() / 2
	final_center_y = cropped_stack.getHeight() / 2
	box_width, box_height = get_rotated_rect_polygon_roi_dims(crop_template)
	IJ.run(cropped_stack, "Specify...", "width=%s height=%s x=%s y=%s centered" % (box_width, box_height, final_center_x, final_center_y))
	cropped_stack_resized = cropped_stack.crop("stack")

	if dataset_metadata_obj.embryo_head_direction == "left":
		IJ.run(cropped_stack_resized, "Rotate 90 Degrees Right", "")
		IJ.run(cropped_stack_resized, "Flip Horizontally", "stack")
	else:
		IJ.run(cropped_stack_resized, "Rotate 90 Degrees Left", "")
	IJ.run(cropped_stack_resized, "Select None", "")
	return cropped_stack_resized

def reset_img_properties(image, voxel_depth):
	"""Reset properties for an ImagePlus image and remove slice lables

	Args:
		image (ImagePlus): input image

	Returns:
		ImagePlus: image with reset properties
	"""
	nslices = image.getNSlices()
	IJ.run(image, "Properties...", "channels=1 slices=%s frames=1 unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=%s.0000 origin=0,0,0" % (nslices, voxel_depth))

	# Equivalent of "Remove Slice Labels" I had to do it manually because
	# "Remove Slice Labels" function always displays the output image
	stack = image.getStack()
	size = image.getStackSize()
	for i in range(1, size + 1):
		stack.setSliceLabel(None, i)
	if size == 1:
		image.setProperty("Label", None)

	return image

def subset_planes(stack_img, planes):
	"""Leave only specified planes in zstack

	Args:
		zstack (ImagePlus): stack to crop
		planes (int, int): (start_plane, end_plane) ends included in the output stack

	Returns:
		ImagePlus: cropped stack
	"""
	# I had to do it manually because the built in function (Keep slices) always displays the image
	stack = stack_img.getStack()
	cropped_stack = ImageStack(stack.getWidth(), stack.getHeight())
	for i in range(planes[0], planes[1] + 1):
		ip = stack.getProcessor(i)
		cropped_stack.addSlice(ip)
	stack_img.setStack(cropped_stack)
	IJ.run(stack_img, "Select None", "")
	return stack_img

def split_montage_image_to_stack(image_processor, nrows, ncols, border_width=0):
	"""Splits a single image with a montage into a stack of images.

	Args:
		image (ImageProcessor): a montage of several images, has to contain a single plane
		nrows (int): number of rows in a montage
		ncols (int): number of columns in a montage
		border_width (int, optional): border width in a montage. Defaults to 0.

	Returns:
		ImageStack: a stack of images from a montage
	"""
	stack = StackMaker.makeStack(StackMaker(), image_processor, nrows, ncols, border_width)
	return stack

def split_montage_stack_to_list_of_stacks(montage_image_stack, nrows, ncols, border_width=0):
	"""Splits a stack of montages into a list of stacks for each image in the montage.

	Args:
		montage_image_stack (ImagePlus): a stack of montages
		nrows (int): number of rows in a montage
		ncols (int): number of columns in a montage
		border_width (int, optional): border width in a montage. Defaults to 0.

	Returns:
		ImagePlus[]: a list of stacks that were used to create montages
	"""
	nstacks = nrows * ncols
	list_of_stacks = [ImageStack() for i in range(nstacks)]
	montage_stack = montage_image_stack.getImageStack()
	for plane_index in range(1, montage_stack.getSize() + 1):
		plane = montage_stack.getProcessor(plane_index)
		stack_from_montage = split_montage_image_to_stack(plane, nrows, ncols, border_width)
		for i, new_stack in enumerate(list_of_stacks):
			new_stack.addSlice(stack_from_montage.getProcessor(i + 1))
	for i, new_stack in enumerate(list_of_stacks):
		list_of_stacks[i] = ImagePlus("New stack from montage %s" % i, new_stack)
	return list_of_stacks

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
		initialX = (stack.getWidth() / 2 - new_width / 2 + extraX) if new_width % 2 == 0 else(stack.getWidth() / 2 -
																							  new_width / 2 + 1 - extraX)
		initialY = (
			stack.getHeight() / 2 - height / 2 + extraY) if height % 2 == 0 else (stack.getHeight() / 2 - height / 2 + 1 - extraY)
		stack.setRoi(initialX, initialY, new_width, height)
		downsampled_imp = stack.crop("stack")
	else:
		IJ.log("You try to upsample the image.  You need an interpolator for that not a downsampler.")
	return downsampled_imp

def make_max_Z_projection_stack_from_folder(input_dir):
	"""Takes all stacks from a directory, does their max-projection, makes a stack of max-projections, saves it to output directory and returns it. 

	Args:
		input_dir (string): path to stacks of images (has to be images from ImagePlus objects)
		output_dir (string): path to output folder

	Returns:
		ImagePlus: stack of max-projections
	"""
	fnames = get_tiffs_in_directory(input_dir)
	if len(fnames) == 0:
		raise Exception("No tiffs to process in %s" % input_dir)
	# Open and stack images
	img_for_dims = IJ.openImage(fnames[0])
	stack_list = []
	for fname in fnames:
		stack = IJ.openImage(fname)
		# Select which dimension to project
		max_proj = z_project_a_stack(stack)
		stack_list.append(max_proj.getProcessor())
	max_proj_stack = ImageStack(img_for_dims.width, img_for_dims.height)
	for slice in stack_list:
		max_proj_stack.addSlice(None, slice)
	max_proj_stack = ImagePlus("max_proj", max_proj_stack)
	max_proj_stack = reset_img_properties(max_proj_stack, voxel_depth=1)
	return max_proj_stack

def get_image_dimensions_from_file(path_to_image):
	image = IJ.openImage(path_to_image)
	return (image.getWidth(), image.getHeight(), image.getStackSize())

def project_image(image_stack, dim_to_project, projection_type):
	data = convertServiceImageUtilsLocal.convert(image_stack, Dataset)

	# Select which dimension to project
	dim = data.dimensionIndex(getattr(Axes, dim_to_project))

	if dim == -1:
		raise Exception("%s dimension not found." % dim_to_project)

	if data.dimension(dim) < 2:
		raise Exception("%s dimension has only one frame." % dim_to_project)

	# Write the output dimensions
	new_dimensions = [data.dimension(d) for d in range(0, data.numDimensions()) if d != dim]

	# Create the output image
	projected = opsServiceImageUtilsLocal.create().img(new_dimensions)

	# Create the op and run it
	proj_op = opsServiceImageUtilsLocal.op(getattr(Ops.Stats, projection_type), data)
	opsServiceImageUtilsLocal.transform().project(projected, data, proj_op, dim)

	projected = opsServiceImageUtilsLocal.run("convert.uint16", projected)

	out = dsServiceImageUtilsLocal.create(projected)
	o = convertServiceImageUtilsLocal.convert(out, ImagePlus)
	ip = o.getProcessor()
	o = ImagePlus("projection", ip)
	return o

def Z_Y_projection_montage_from_image_stack(stack):
	z_proj = project_image(stack, "Z", "Max")
	y_proj = project_image(stack, "Y", "Max")
	combined = StackCombiner.combineVertically(StackCombiner(), z_proj.getStack(), y_proj.getStack())
	return ImagePlus("montage", combined)


###### Image histogram/intensity management

def get_histogram_thresholds(image_processor, percent_overexposed_pixels=1):
	"""Get upper and lower thresholds for the histogram of the image for contrast adjustment. 
	Lower histogram is determined by Triangle thresholding algorithm.
	Upper histogram is determined in such a way that globally determined persentage of pixels in the image after applying lower threshold will be overexposed.

	Args:
		image_processor (ImageProcessor): image to determine histogram thresholds from
        percent_overexposed_pixels (int): percentage of pixels that will be overexposed in the adjusted image


	Returns:
		(int, int): a tuple of lower and upper thresholds
	"""
	hist = image_processor.getHistogram()
	lower_threshold = Auto_Threshold.Triangle(hist)
	num_overexposed_pixels = 0
	upper_threshold = 65535
	sum_elem_above_threshold = sum(hist[i] for i in range(lower_threshold, len(hist)))

	# making it so 1% of all pixels will be overexposed
	num_overexposed_pixels_threshold = float(percent_overexposed_pixels) / 100 * sum_elem_above_threshold
	for i, num_pixels_with_this_value in reversed(list(enumerate(hist))):
		num_overexposed_pixels += num_pixels_with_this_value
		if num_overexposed_pixels > num_overexposed_pixels_threshold:
			upper_threshold = i
			break
	return (lower_threshold, upper_threshold)

def auto_contrast_by_histogram_thresholds(image, thresholds):
	"""Adjust contrast of the image by removing pixels based on the lower and upper intensity thresholds.

	Args:
		image (ImagePlus): image to adjust
		thresholds (int, int): a tuple of lower and upper thresholds
	"""
	ip = image.getProcessor()
	lower_threshold, upper_threshold = thresholds

	# logging.info("Chosen theese values to adjust image histogram, min: %s max: %s" % (lower_threshold, upper_threshold))
	image.setDisplayRange(lower_threshold, upper_threshold)
	IJ.run(image, "Apply LUT", "")

def threshold_histogram_stack(imp_stack, percent_overexposed_pixels=1):
	"""Adjust contrast of the image stack by removing pixels based on the lower and upper intensity thresholds. Intensity thresholds are determined from the middle plane in the stack by get_histogram_thresholds function. 

	Args:
		imp_stack (ImagePlus): image stack to adjust
        percent_overexposed_pixels (int): percentage of pixels that will be overexposed in the adjusted image

	Raises:
		Exception: Not a 16-bit image has been provided.

	Returns:
		(ImagePlus, (int, int)): a tuple of adjusted stack of images and a tuple of lower and upper histogram thresholds used for adjustment
	"""
	if imp_stack.getBitDepth() != 16:
		raise Exception("Only 16-bit images can be auto histogram adjusted.")
	adjusted_stack = []
	stack = imp_stack.getStack()
	nplanes = imp_stack.getNSlices()
	middle_plane = int(round(float(nplanes) / 2))
	histogram_thresholds = get_histogram_thresholds(stack.getProcessor(middle_plane), percent_overexposed_pixels)

	for i in range(1, nplanes + 1):
		imp2 = ImagePlus("slice_iterator", stack.getProcessor(i))
		auto_contrast_by_histogram_thresholds(imp2, histogram_thresholds)
		adjusted_stack.append(imp2)
	adjusted_stack = ImagePlus(
		"Adjusted contrast", ImageStack.create(adjusted_stack))
	return (adjusted_stack, histogram_thresholds)

def match_histograms_stack(stack):
	"""Adjust histrograms for the provided image stack. Mutates provided stack. Resulting histograms for each plane have the same area
	under the curve and look equally bright.

	Args:
		stack (ImagePlus): a stack of images
	"""
	bleach_obj = BleachCorrection_MH(stack)
	bleach_obj.doCorrection()

def get_median_intensity(image):
	image.resetRoi()
	IJ.run("Clear Results")
	table = ResultsTable()
	analyzer = Analyzer(image, Measurements.MEDIAN, table)
	analyzer.measure()
	median = round(table.getValue("Median", 0), ndigits=0)
	return median

def get_mean_intensity(image):
	image.resetRoi()
	IJ.run("Clear Results")
	table = ResultsTable()
	analyzer = Analyzer(image, Measurements.MEAN, table)
	analyzer.measure()
	mean = round(table.getValue("Mean", 0), ndigits=0)
	return mean

def get_zero_pixels_mask(image_16bit):
	ip = image_16bit.getProcessor().duplicate()  # get pixel array?, as a copy
	mask = ImagePlus("mask", ip)
	hist = ip.getHistogram()
	ip.setThreshold(0, 1, ImageProcessor.NO_LUT_UPDATE)
	IJ.run(mask, "Convert to Mask", "")
	return mask

def fill_zero_pixels_with_median_intensity(image):
	# This check has to be here otherwize the image will get flooded with median value, and thresholds do not work
	if image.getProcessor().getStats().min > 0:
		return
	median = get_median_intensity(image)
	zero_mask = get_zero_pixels_mask(image)
	IJ.run(zero_mask, "Create Selection", "")
	image.setRoi(zero_mask.getRoi())
	# when selection is inverted you do not cover a line of pixels on selection border, so it has to be expanded.
	RoiEnlarger.enlarge(image, 2);
	IJ.run(image, "Set...", "value=%s" % median)
	image.resetRoi()


###### Roi managing

def polygon_to_rotated_rect_roi(roi):
	if isinstance(roi, RotatedRectRoi):
		return roi
	if roi.getNCoordinates() != 4:
		raise Exception(
			"polygon_to_rotated_rect_roi Can only convert rectangles. Received polygon with %s points." % roi.getNCoordinates())
	x1 = roi.getFloatPolygon().xpoints[0]
	x2 = roi.getFloatPolygon().xpoints[1]
	x3 = roi.getFloatPolygon().xpoints[2]
	x4 = roi.getFloatPolygon().xpoints[3]
	y1 = roi.getFloatPolygon().ypoints[0]
	y2 = roi.getFloatPolygon().ypoints[1]
	y3 = roi.getFloatPolygon().ypoints[2]
	y4 = roi.getFloatPolygon().ypoints[3]
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

def get_rectangular_polygon_roi_angle(roi):
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

def get_rotated_rect_polygon_roi_dims(roi):
	x1 = roi.getFloatPolygon().xpoints[0]
	x2 = roi.getFloatPolygon().xpoints[1]
	x3 = roi.getFloatPolygon().xpoints[2]
	y1 = roi.getFloatPolygon().ypoints[0]
	y2 = roi.getFloatPolygon().ypoints[1]
	y3 = roi.getFloatPolygon().ypoints[2]
	dim1 = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
	dim2 = math.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)
	width = int(round(max(dim1, dim2)))
	height = int(round(min(dim1, dim2)))
	# width -= width % 4 
	# height -= height % 4 
	# logging.info("Determined rectangular roi dims before rounding: %s, %s" % (max(dim1, dim2), min(dim1, dim2)))
	# logging.info("\tRectangular roi dims after rounding: %s, %s" % (width, height))
	return width, height

def save_roi(roi, roi_save_path):
	roim = RoiManager(False)
	roim.runCommand('reset')
	roim.addRoi(roi)
	roim.select(0)
	roim.runCommand("Save", roi_save_path)
	roim.close()

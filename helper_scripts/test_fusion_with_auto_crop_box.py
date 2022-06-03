# @ File(label='Input directory with all datasets', style='directory') input_datasets_dir
# @ File(label='Results directory', style='directory') results_dir
# @ Boolean (label='Use previously fused datasets?', value=false) use_fused_cache
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert

# Script downsamples first timepoint from each direction from the raw images folder. Fuses the downsampled images, makes Z and Y projections calculates embryo bounding box from these projections and fuses the downsampled 1 timepoint again with new crop box.

import os
import re
import math
from distutils.dir_util import mkpath
import shutil
import xml.etree.ElementTree as ET


from java.io import File
from ij.io import FileSaver
from ij import IJ, ImagePlus, ImageStack, WindowManager
from net.imagej.axis import Axes
from net.imagej.ops import Ops
from net.imglib2.view import Views
from net.imagej import Dataset
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
from ij.plugin.filter import Analyzer
from ij.plugin import RoiEnlarger

RAW_IMAGES_DIR_NAME = "(P0)-ZStacks-Raw"

class FredericFile:
	"""
	File naming for light-sheet image files. Frederic Strobl™ 
	Example file name: MGolden2022A-DS0001TP0001DR0001CH0001PL(ZS).tif
	:param file_name: full image file name
	:type file_name: str
	"""
	dataset_name = ""
	dataset_id = 0
	time_point = ""
	direction = ""
	channel = ""
	plane = ""
	extension = ""

	def __init__(self, file_name):
		name_parts = re.split("-DS|TP|DR|CH|PL|\.", file_name)
		if len(name_parts) != 7:
			raise Exception(
				"Image file name is improperly formatted! Check documentation inside the script.")
		self.dataset_name, self.dataset_id, self.time_point, self.direction, self.channel, self.plane, self.extension = name_parts
		self.extension.lower()

	def get_name(self):
		return "%s-DS%sTP%sDR%sCH%sPL%s.%s" % (self.dataset_name,
                                         self.dataset_id,
                                         self.time_point,
                                         self.direction,
                                         self.channel,
                                         self.plane,
                                         self.extension)


def get_tiffs_in_directory(directory):
	"""Get all TIFF file paths in a directory. Subdirectories are not searched.

	Args:
		directory (str): full path to directory

	Raises:
		Exception: If no images were found

	Returns:
		str[]: list of full paths to tiff files
	"""
	file_names = []
	for fname in os.listdir(directory):
		if fname.lower().endswith(".tif"):
			file_names.append(os.path.join(directory, fname))
	file_names = sorted(file_names)
	return file_names

def get_image_dimensions_from_file(path_to_image):
    image = IJ.openImage(path_to_image)
    return (image.getWidth(), image.getHeight(), image.getStackSize())

def create_and_fuse_dataset(dataset_dir, file_pattern, timepoints, angles_from_to, view_stack_dims):
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

    IJ.run(
        "Define Bounding Box for Fusion",
        "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] bounding_box=[Maximal Bounding Box spanning all transformed views] bounding_box=[Maximal Bounding Box spanning all transformed views] bounding_box_name=max_box minimal_x=0 minimal_y=0 minimal_z=0 maximal_x=%s maximal_y=%s maximal_z=%s"
        % (dataset_xml, view_stack_dims[0], view_stack_dims[1], view_stack_dims[2]))

    print("saving fused dataset to: ", output_fused_path)
    IJ.run(
        "Fuse dataset ...",
        "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] bounding_box=max_box downsampling=1 pixel_type=[16-bit unsigned integer] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend produce=[Each timepoint & channel] fused_image=[Save as new XML Project (TIFF)] export_path=%s" % (dataset_xml, output_fused_path))


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


def project_image(image_stack, dim_to_project, projection_type):
    data = convert.convert(image_stack, Dataset)

    # Select which dimension to project
    dim = data.dimensionIndex(getattr(Axes, dim_to_project))

    if dim == -1:
        raise Exception("%s dimension not found." % dim_to_project)

    if data.dimension(dim) < 2:
        raise Exception("%s dimension has only one frame." % dim_to_project)

    # Write the output dimensions
    new_dimensions = [data.dimension(d) for d in range(0, data.numDimensions()) if d != dim]

    # Create the output image
    projected = ops.create().img(new_dimensions)

    # Create the op and run it
    proj_op = ops.op(getattr(Ops.Stats, projection_type), data)
    ops.transform().project(projected, data, proj_op, dim)

    projected = ops.run("convert.uint16", projected)

    out = ds.create(projected)
    o = convert.convert(out, ImagePlus)
    ip = o.getProcessor()
    o = ImagePlus("projection", ip)
    return o


def make_threshold_montage(image):
    IJ.run(image, "Auto Threshold", "method=[Try all] white")
    montage = WindowManager.getImage("Montage")
    return montage


# Functions for segmenting embryo, cropping and fusing cropped

def get_median_intensity(image):
    image.resetRoi()
    IJ.run("Clear Results")
    table = ResultsTable()
    analyzer = Analyzer(image, Measurements.MEDIAN, table)
    analyzer.measure()
    median = round(table.getValue("Median", 0), ndigits=0)
    print("Median value is %s" % median)
    return median


def get_mean_intensity(image):
    image.resetRoi()
    IJ.run("Clear Results")
    table = ResultsTable()
    analyzer = Analyzer(image, Measurements.MEAN, table)
    analyzer.measure()
    mean = round(table.getValue("Mean", 0), ndigits=0)
    print("Mean value is %s" % mean)
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


def get_embryo_bounding_rectangle_and_params(max_projection, show_mask=False, fill_zeros=True):
    if fill_zeros:
        fill_zero_pixels_with_median_intensity(max_projection)

    # image_hist = max_projection.getProcessor().getHistogram()
    # mean_intensity = get_mean_intensity(max_projection)
    # median_intensity = get_median_intensity(max_projection)

    mask = get_embryo_mask(image_16bit=max_projection, threshold_type="triangle")
    # threshold_choosing_intensity =  int(1.4 * median_intensity)
    # if mean_intensity < threshold_choosing_intensity:
    #     print("Image average intesity < %s using triangle threshold." % threshold_choosing_intensity)
    #     mask = get_embryo_mask(image_16bit=max_projection, threshold_type="triangle")
    # else:
    #     print("Image average intesity > %s using triangle threshold." % threshold_choosing_intensity)
    #     mask = get_embryo_mask(image_16bit=max_projection, threshold_type="triangle")   
    # else:
    #     print("Image average intesity > %s using huang2 threshold." % threshold_choosing_intensity)
    #     mask = get_embryo_mask(image_16bit=max_projection, threshold_type="huang2")   

    if show_mask:
        mask.show() 

    num_mask_smoothening_iterations = int(round(mask.getWidth() * 0.028))
    # removing dirt, fixing edges
    smoothen_embryo_mask(mask, num_mask_smoothening_iterations)

    bounding_roi_stuff = bounding_roi_and_params_from_embryo_mask(mask=mask, spacing_around_embryo=8)
    bounding_roi_stuff["mask"] = mask
    return bounding_roi_stuff


def project_image(image, dim_to_project, projection_type):
    data = convert.convert(image, Dataset)


    # Select which dimension to project
    dim = data.dimensionIndex(getattr(Axes, dim_to_project))

    if dim == -1:
        raise Exception("%s dimension not found." % dim_to_project)

    if data.dimension(dim) < 2:
        raise Exception("%s dimension has only one frame." % dim_to_project)

    # Write the output dimensions
    new_dimensions = [data.dimension(d) for d in range(0, data.numDimensions()) if d != dim]

    # Create the output image
    projected = ops.create().img(new_dimensions)


    # Create the op and run it
    proj_op = ops.op(getattr(Ops.Stats, projection_type), data)
    ops.transform().project(projected, data, proj_op, dim)

    projected = ops.run("convert.uint16", projected)

    out = ds.create(projected)
    o = convert.convert(out, ImagePlus)
    ip = o.getProcessor()
    o = ImagePlus("projection", ip)
    return o


def get_fusion_tranformation_from_xml_file(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    for transform in root.iter('ViewRegistrations'):
        # This is very dependent on the BigStitcher XML fromat being consistent and not changing! 
        if transform[0][0][0].text == "fusion bounding box":
            matrix_str = transform[0][0].find("affine").text.split(" ")
            
            matrix = [[0 for col in range(4)] for row in range(4)]
            matrix[3][3] = 1

            # iterate through rows
            for i in range(3):
                # iterate through columns
                for j in range(4):
                    matrix[i][j] = float(matrix_str[i * 4 + j])          
            return matrix
    return False


def rotate_bigstitcher_dataset(dataset_xml_path, axis_of_rotation, angle, timepoint=0):
    print("Rotating dataset around: %s for angle=%s" % (axis_of_rotation, angle))
    IJ.run("Apply Transformations", "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Rigid apply=[Current view transformations (appends to current transforms)] define=[Rotation around axis] same_transformation_for_all_angles axis_timepoint_%s_channel_0_illumination_0_all_angles=%s-axis rotation_timepoint_%s_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, timepoint,  axis_of_rotation, timepoint, angle))


def apply_transformation_bigstitcher_dataset(dataset_xml_path, affine_matrix, timepoint=0):
    short_affine_matrix = [0 for i in range(12)]
    for i in range(len(affine_matrix) - 1):
        for j in range(len(affine_matrix[0])):
            short_affine_matrix[i * 4 + j] = affine_matrix[i][j]
    print("Applying transformation with matrix %s" % short_affine_matrix)

    IJ.run("Apply Transformations", "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Affine apply=[Current view transformations (appends to current transforms)] same_transformation_for_all_angles timepoint_%s_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, timepoint, short_affine_matrix))


def define_bounding_box_for_fusion(dataset_xml_path, box_coords, bounding_box_name, timepoint=0):
    IJ.run("Define Bounding Box for Fusion", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] processing_timepoint=[Timepoint %s] bounding_box=[Maximal Bounding Box spanning all transformed views] bounding_box_name=%s minimal_x=%s minimal_y=%s minimal_z=%s maximal_x=%s maximal_y=%s maximal_z=%s" % (dataset_xml_path, timepoint, bounding_box_name, box_coords["x_min"], box_coords["y_min"], box_coords["z_min"], box_coords["x_max"], box_coords["y_max"], box_coords["z_max"]))


def fuse_dataset_to_display(dataset_xml_path, bounding_box_name):
    IJ.run("Fuse dataset ...", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] bounding_box=%s downsampling=1 pixel_type=[16-bit unsigned integer] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend produce=[Each timepoint & channel] fused_image=[Display using ImageJ]" % (dataset_xml_path, bounding_box_name))


def fuse_dataset_to_tiff(dataset_xml_path, bounding_box_name, fused_xml_path):
    IJ.run(
    "Fuse dataset ...",
    "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] bounding_box=%s downsampling=1 pixel_type=[16-bit unsigned integer] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend produce=[Each timepoint & channel] fused_image=[Save as new XML Project (TIFF)] export_path=%s" % (dataset_xml_path, bounding_box_name, fused_xml_path))


def get_4_4_identity_matrix():
    identity_matrix = [[0 for col in range(4)] for row in range(4)]
    identity_matrix[0][0] = 1
    identity_matrix[1][1] = 1
    identity_matrix[2][2] = 1
    identity_matrix[3][3] = 1
    return identity_matrix


def multiply_matrices(X, Y):
	C = [[0 for col in range(len(Y[0]))] for row in range(len(X))]
	print(C)
	# iterate through rows of X
	for i in range(len(X)):
		# iterate through columns of Y
		for j in range(len(Y[0])):
			# iterate through rows of Y
			for k in range(len(Y)):
				C[i][j] += X[i][k] * Y[k][j]

	return C


def inverse_matrix(m):
    def return_transpose(mat):
        return map(list,zip(*mat))

    def return_matrix_minor(mat,i,j):
        return [row[:j] + row[j+1:] for row in (mat[:i]+mat[i+1:])]

    def return_determinant(m):
        if len(m) == 2:
            return m[0][0]*m[1][1]-m[0][1]*m[1][0]

        determinant = 0
        for c in range(len(m)):
            determinant += ((-1)**c)*m[0][c]*return_determinant(return_matrix_minor(m,0,c))
        return determinant
    
    determinant = return_determinant(m)
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    cfs = []
    for r in range(len(m)):
        cfRow = []
        for c in range(len(m)):
            minor = return_matrix_minor(m,r,c)
            cfRow.append(((-1)**(r+c)) * return_determinant(minor))
        cfs.append(cfRow)
    cfs = return_transpose(cfs)
    for r in range(len(cfs)):
        for c in range(len(cfs)):
            cfs[r][c] = cfs[r][c]/determinant
    return cfs


def save_roi(roi, roi_save_path):
    roim = RoiManager(False)
    roim.runCommand('reset')
    roim.addRoi(roi)
    roim.select(0)
    roim.runCommand("Save", roi_save_path)
    roim.close()


def save_tiff(image, path, compress_on_save=False):
	if os.path.exists(path):
		os.remove(path)
	if compress_on_save == False:
		fs = FileSaver(image)
		fs.saveAsTiff(path)
	else:
		IJ.run(image, "Bio-Formats Exporter", "save=%s export compression=zlib" % path)


def Z_Y_projection_montage_from_image_stack(stack):
    z_proj = project_image(stack, "Z", "Max")
    y_proj = project_image(stack, "Y", "Max")
    combined = StackCombiner.combineVertically(StackCombiner(), z_proj.getStack(), y_proj.getStack())
    return ImagePlus("montage", combined)


def segment_embryo_and_fuse_again_cropping_around_embryo(raw_dataset_xml_path, fused_xml_path, z_projection, y_projection):
    dataset_dir = os.path.dirname(raw_dataset_xml_path)
    cropped_fused_dir = os.path.join(dataset_dir, "cropped_fusion")
    if not os.path.exists(cropped_fused_dir):
        mkpath(cropped_fused_dir)

    z_bounding_roi = get_embryo_bounding_rectangle_and_params(max_projection=z_projection, show_mask=False)
    y_bounding_roi = get_embryo_bounding_rectangle_and_params(max_projection=y_projection, show_mask=False)

    if abs(z_bounding_roi["embryo_center_x"] - y_bounding_roi["embryo_center_x"]) > 5:
        print("Could not reliably determine embryo X center position from Z and Y projections.")

    print("Embryo box params from Z projection:\n%s" % z_bounding_roi)
    print("Embryo box params from Y projection:\n%s " % y_bounding_roi)

    
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
        exit(1)
    print("transft to raw %s" % transf_to_raw)
    transformation_embryo_to_center = multiply_matrices(transformation_embryo_to_center, transf_to_raw)
    transformation_embryo_to_center = inverse_matrix(transformation_embryo_to_center)
    print("transformation_embryo_to_center matrix %s" % transformation_embryo_to_center)

    apply_transformation_bigstitcher_dataset(raw_dataset_xml_path, transformation_embryo_to_center)
    embryo_crop_box = {
        "x_min" : -1 * int(z_bounding_roi["embryo_length"] / 2 + 5),
        "y_min" : -1 * int(z_bounding_roi["embryo_width"] / 2 + 5),
        "z_min" : -1 * int(y_bounding_roi["embryo_width"] / 2 + 5),
        "x_max" : int(z_bounding_roi["embryo_length"] / 2 + 5),
        "y_max" : int(z_bounding_roi["embryo_width"] / 2 + 5),
        "z_max" : int(y_bounding_roi["embryo_width"] / 2 + 5)
    }

    rotate_bigstitcher_dataset(raw_dataset_xml_path, "z", round(z_bounding_roi["bounding_rect_angle"], 1))
    rotate_bigstitcher_dataset(raw_dataset_xml_path, "y", -1 * round(y_bounding_roi["bounding_rect_angle"], 1))


    define_bounding_box_for_fusion(raw_dataset_xml_path, embryo_crop_box, "embryo_cropped")
    cropped_fused_path = os.path.join(cropped_fused_dir, "cropped_fusion.xml")
    cropped_fusion_stack_path = os.path.join(cropped_fused_dir, "fused_tp_0_ch_0.tif")
    fuse_dataset_to_tiff(raw_dataset_xml_path, "embryo_cropped", cropped_fused_path)
    zy_cropped_projections = Z_Y_projection_montage_from_image_stack(IJ.openImage(cropped_fusion_stack_path))

    return zy_cropped_projections

# Main process

def process_dataset(input_dataset_dir, results_dir):

    raw_img_dir = os.path.join(input_dataset_dir, RAW_IMAGES_DIR_NAME, "CH0001", "DR0001")
    if not os.path.exists(raw_img_dir):
        return

    tiffs = get_tiffs_in_directory(raw_img_dir)
    if tiffs == []:
        return
    last_file = os.path.basename(tiffs[-1])
    print(last_file)
    last_timepoint_name = FredericFile(last_file)
    timepoints_to_downsample = [1]
    dataset_name = last_timepoint_name.dataset_name + "_" + os.path.basename(input_dataset_dir)

    direction_dirs = [os.path.join(input_dataset_dir, RAW_IMAGES_DIR_NAME, "CH0001", "DR000%s" % i) for i in range(1, 5)]

    results_subdir = os.path.join(results_dir, dataset_name)
    if not os.path.exists(results_subdir):
        mkpath(results_subdir)
    
    cropped_projections_subdir = os.path.join(results_dir, "cropped_fused_projections")
    if not os.path.exists(cropped_projections_subdir):
        mkpath(cropped_projections_subdir)

    view_image_name = last_timepoint_name

    tp = 1
    for direction, direction_dir in zip(range(1, 5), direction_dirs):
        view_image_name.time_point = "%04d" % tp
        view_image_name.direction = "%04d" % direction
        downsampled_path = os.path.join(results_dir, dataset_name, view_image_name.get_name())
        if os.path.exists(downsampled_path):
            print("Found existing file %s, skipping downsampling." % downsampled_path)
            continue

        full_image_path = os.path.join(direction_dir, view_image_name.get_name())
        full_image = IJ.openImage(full_image_path)
        downsampled_image = downsample_image_stack(full_image, full_image.getWidth() / 4)
        print("Saving downsampled image to %s" % downsampled_path)
        fs = FileSaver(downsampled_image)
        fs.saveAsTiff(os.path.join(results_dir, dataset_name, view_image_name.get_name()))
    fuse_file_pattern = last_timepoint_name
    fuse_file_pattern.time_point = "0001"
    view_dims = get_image_dimensions_from_file(os.path.join(results_subdir, fuse_file_pattern.get_name()))
    fuse_file_pattern.direction = "{aaaa}"
    if os.path.exists(os.path.join(results_subdir, "fused.xml")) and use_fused_cache == True:
        print("Skipping fusion for %s, found exiting fused.xml file." % results_subdir)
    else:
        print("Fusing: " + fuse_file_pattern.get_name())
        create_and_fuse_dataset(results_subdir, fuse_file_pattern.get_name(), timepoints_to_downsample, (1, 4), view_dims)

    fused_file_path = os.path.join(results_subdir, "fused_tp_0_ch_0.tif")

    fused_stack = IJ.openImage(fused_file_path)
    z_projection = project_image(fused_stack, "Z", "Max")
    y_projection = project_image(fused_stack, "Y", "Max")
    # RankFilters().rank(z_projection.getProcessor(), 4, RankFilters.MEDIAN)
    # RankFilters().rank(y_projection.getProcessor(), 4, RankFilters.MEDIAN)

    # fs = FileSaver(z_projection)
    # fs.saveAsTiff(os.path.join(results_subdir, "tp_%s_z_max_projection.tiff" % tp))
    # fs = FileSaver(y_projection)
    # fs.saveAsTiff(os.path.join(results_subdir, "tp_%s_y_max_projection.tiff" % tp))

    # z_montage = make_threshold_montage(z_projection)
    # z_montage_resized = z_montage.resize(1800, 1200, "bilinear")

    # fs = FileSaver(z_montage_resized)
    # output_z_montage_name = view_image_name
    # output_z_montage_name.dataset_name = dataset_name + "_Z_threshold_montage"
    # output_z_montage_name.direction = "(ZM)"
    # output_z_montage_name.time_point = "%04d" % tp
    # fs.saveAsTiff(os.path.join(results_dir, output_z_montage_name.get_name()))
    # print("Saving threshold montage to: %s" % os.path.join(results_dir, output_z_montage_name.get_name()))
    # z_montage.close()

    # y_montage = make_threshold_montage(y_projection)
    # y_montage_resized = y_montage.resize(1800, 1200, "bilinear")
    # fs = FileSaver(y_montage_resized)
    # output_y_montage_name = view_image_name
    # output_y_montage_name.dataset_name = dataset_name + "_Y_threshold_montage"
    # output_y_montage_name.direction = "(ZM)"
    # output_y_montage_name.time_point = "%04d" % tp
    # fs.saveAsTiff(os.path.join(results_dir, output_y_montage_name.get_name()))
    # y_montage.close()
    shutil.copy2(os.path.join(results_subdir, "dataset.xml"), os.path.join(results_subdir, "dataset_centered.xml"))
    zy_cropped_projections = segment_embryo_and_fuse_again_cropping_around_embryo(
                                                        raw_dataset_xml_path=os.path.join(results_subdir, "dataset_centered.xml"),
                                                        fused_xml_path=os.path.join(results_subdir, "fused.xml"),
                                                        z_projection=z_projection,
                                                        y_projection=y_projection)
    IJ.setBackgroundColor(255, 255, 255)
    IJ.run(zy_cropped_projections, "Canvas Size...", "width=350 height=350 position=Center")                                            
    save_tiff(zy_cropped_projections, os.path.join(cropped_projections_subdir, dataset_name + "_zy_cropped_fused_projections.tiff"))

input_datasets_dir = input_datasets_dir.getAbsolutePath()
results_dir = results_dir.getAbsolutePath()

dirs = os.listdir(input_datasets_dir)
for dataset_dir in dirs:
    dataset_dir = os.path.join(input_datasets_dir, dataset_dir)
    print("dataset_dir", dataset_dir)
    if os.path.isdir(dataset_dir):
        process_dataset(dataset_dir, results_dir)

print("DONE")

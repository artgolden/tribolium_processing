# @ File(label='Input directory with all datasets', style='directory') input_datasets_dir
# @ File(label='Results directory', style='directory') results_dir
# @ Boolean (label='Use previously fused datasets?', value=false) use_fused_cache
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert

# Script downsamples first, middle and last timepoints from each direction from the raw images folder. Fuses the downsampled images, makes Z and Y projections for each timepoint and produces a threshold montage for each projection.

import os
import re
import math
from distutils.dir_util import mkpath
import xml.etree.ElementTree as ET


from java.io import File
from ij.io import FileSaver
from ij import IJ, ImagePlus, ImageStack, WindowManager
from net.imagej.axis import Axes
from net.imagej.ops import Ops
from net.imglib2.view import Views
from net.imagej import Dataset
from ij.plugin.filter import RankFilters



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


def get_bounding_box_coords_from_file(xml_path, box_name):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    for box in root.iter('BoundingBoxDefinition'):
         if box.attrib["name"] == box_name:
             min = box.find("min").text.split(" ")
             max = box.find("max").text.split(" ")
             box_coords = {
                 "x_min" : int(min[0]),
                 "y_min" : int(min[1]),
                 "z_min" : int(min[2]),
                 "x_max" : int(max[0]),
                 "y_max" : int(max[1]),
                 "z_max" : int(max[2])
             }
    return box_coords


def fuse_dataset(dataset_dir, file_pattern, timepoints, angles_from_to):
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
        "define_dataset=[Manual Loader (TIFF only, ImageJ Opener)] project_filename=%s multiple_timepoints=[YES (one file per time-point)] multiple_channels=[NO (one channel)] _____multiple_illumination_directions=[NO (one illumination direction)] multiple_angles=[YES (one file per angle)] multiple_tiles=[NO (one tile)] image_file_directory=%s image_file_pattern=%s timepoints_=%s acquisition_angles_=%s-%s calibration_type=[Same voxel-size for all views] calibration_definition=[User define voxel-size(s)] imglib2_data_container=[ArrayImg (faster)] pixel_distance_x=1.00000 pixel_distance_y=1.00000 pixel_distance_z=1.00000 pixel_unit=um" % (dataset_xml_name, dataset_dir, file_pattern, timepoints, angles_from_to[0], angles_from_to[1]))

    IJ.run(
        "Detect Interest Points for Pairwise Registration",
        "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] type_of_interest_point_detection=Difference-of-Gaussian label_interest_points=beads subpixel_localization=[3-dimensional quadratic fit] interest_point_specification=[Advanced ...] downsample_xy=[Match Z Resolution (less downsampling)] downsample_z=1x use_same_min sigma=%s threshold=%s find_maxima maximum_number=%s type_of_detections_to_use=Brightest compute_on=[CPU (Java)]"
        % (dataset_xml, sigma, threshold, max_detections_per_view))

    IJ.run(
        "Register Dataset based on Interest Points",
        "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] registration_algorithm=[Fast descriptor-based (rotation invariant)] registration_over_time=[Register timepoints individually] registration_in_between_views=[Only compare overlapping views (according to current transformations)] interest_points=beads fix_views=[Fix first view] map_back_views=[Do not map back (use this if views are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 redundancy=%s significance=%s allowed_error_for_ransac=%s number_of_ransac_iterations=%s"
        % (dataset_xml, redundancy, significance, allowed_error_for_ransac, number_of_ransac_iterations))

    IJ.run(
        "Register Dataset based on Interest Points",
        "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] registration_algorithm=[Fast descriptor-based (rotation invariant)] registration_over_time=[Match against one reference timepoint (no global optimization)] registration_in_between_views=[Only compare overlapping views (according to current transformations)] interest_points=beads reference=%s consider_each_timepoint_as_rigid_unit transformation=Translation regularize_model model_to_regularize_with=Rigid lamba=0.10 redundancy=1 significance=10 allowed_error_for_ransac=2 number_of_ransac_iterations=Normal interestpoint_grouping=[Group interest points (simply combine all in one virtual view)] interest=5" % (dataset_xml, reference_timepoint))


    IJ.run(
        "Define Bounding Box for Fusion",
        "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] bounding_box=[Maximal Bounding Box spanning all transformed views] bounding_box=[Maximal Bounding Box spanning all transformed views] bounding_box_name=max_box"
        % (dataset_xml))

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


def process_dataset(input_dataset_dir, results_dir):

    raw_img_dir = os.path.join(input_dataset_dir, RAW_IMAGES_DIR_NAME, "CH0001", "DR0001")

    tiffs = get_tiffs_in_directory(raw_img_dir)
    last_file = os.path.basename(tiffs[-1])
    print(last_file)
    last_timepoint_name = FredericFile(last_file)
    ntimepoints = int(last_timepoint_name.time_point)
    timepoints_to_downsample = [1, ntimepoints / 2, ntimepoints]
    dataset_name = last_timepoint_name.dataset_name + "_" + os.path.basename(input_dataset_dir)

    direction_dirs = [os.path.join(input_dataset_dir, RAW_IMAGES_DIR_NAME, "CH0001", "DR000%s" % i) for i in range(1, 5)]

    results_subdir = os.path.join(results_dir, dataset_name)
    if not os.path.exists(results_subdir):
        mkpath(results_subdir)

    view_image_name = last_timepoint_name

    for tp in timepoints_to_downsample:
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
    fuse_file_pattern.direction = "{aaaa}"
    fuse_file_pattern.time_point = "{tttt}"
    if os.path.exists(os.path.join(results_subdir, "fused.xml")) and use_fused_cache == True:
        print("Skipping fusion for %s, found exiting fused.xml file." % results_subdir)
    else:
        print("Fusing: " + fuse_file_pattern.get_name())
        fuse_dataset(results_subdir, fuse_file_pattern.get_name(), timepoints_to_downsample, (1, 4))

    for tp in timepoints_to_downsample:
        fused_file_path = os.path.join(results_subdir, "fused_tp_%s_ch_0.tif" % tp)

        fused_stack = IJ.openImage(fused_file_path)
        z_projection = project_image(fused_stack, "Z", "Max")
        y_projection = project_image(fused_stack, "Y", "Max")
        RankFilters().rank(z_projection.getProcessor(), 4, RankFilters.MEDIAN)
        RankFilters().rank(y_projection.getProcessor(), 4, RankFilters.MEDIAN)

        fs = FileSaver(z_projection)
        fs.saveAsTiff(os.path.join(results_subdir, "tp_%s_z_max_projection.tiff" % tp))
        fs = FileSaver(y_projection)
        fs.saveAsTiff(os.path.join(results_subdir, "tp_%s_y_max_projection.tiff" % tp))

        z_montage = make_threshold_montage(z_projection)
        z_montage_resized = z_montage.resize(1800, 1200, "bilinear")

        fs = FileSaver(z_montage_resized)
        output_z_montage_name = view_image_name
        output_z_montage_name.dataset_name = dataset_name + "_Z_threshold_montage"
        output_z_montage_name.direction = "(ZM)"
        output_z_montage_name.time_point = "%04d" % tp
        fs.saveAsTiff(os.path.join(results_dir, output_z_montage_name.get_name()))
        print("Saving threshold montage to: %s" % os.path.join(results_dir, output_z_montage_name.get_name()))
        z_montage.close()

        y_montage = make_threshold_montage(y_projection)
        y_montage_resized = y_montage.resize(1800, 1200, "bilinear")
        fs = FileSaver(y_montage_resized)
        output_y_montage_name = view_image_name
        output_y_montage_name.dataset_name = dataset_name + "_Y_threshold_montage"
        output_y_montage_name.direction = "(ZM)"
        output_y_montage_name.time_point = "%04d" % tp
        fs.saveAsTiff(os.path.join(results_dir, output_y_montage_name.get_name()))
        y_montage.close()

input_datasets_dir = input_datasets_dir.getAbsolutePath()
results_dir = results_dir.getAbsolutePath()

dirs = os.listdir(input_datasets_dir)
for dataset_dir in dirs:
    dataset_dir = os.path.join(input_datasets_dir, dataset_dir)
    print("dataset_dir", dataset_dir)
    if os.path.isdir(dataset_dir):
        process_dataset(dataset_dir, results_dir)

print("DONE")

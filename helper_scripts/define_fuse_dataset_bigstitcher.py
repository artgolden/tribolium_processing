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

    max_box_coords = get_bounding_box_coords_from_file(dataset_xml, "max_box")


    print("saving fused dataset to: ", output_fused_path)
    IJ.run(
        "Fuse dataset ...",
        "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] bounding_box=max_box downsampling=1 pixel_type=[16-bit unsigned integer] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend produce=[Each timepoint & channel] fused_image=[Save as new XML Project (TIFF)] export_path=%s" % (dataset_xml, output_fused_path))

fuse_dataset("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_dataset", "MGolden2022A-DS0016TP{tttt}DR{aaaa}CH0001PL(ZS).tif", [1], (1, 4))


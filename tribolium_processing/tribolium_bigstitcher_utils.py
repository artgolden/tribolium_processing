import logging
import os
import sys
import time
import xml.etree.ElementTree as ET



from ij import IJ, WindowManager
from ij.io import FileSaver

from java.lang.System import getProperty
sys.path.append(os.path.join(getProperty("fiji.dir"), "plugins", "tribolium_processing")) 
from tribolium_file_utils import *


# small Utils

def save_tiff_simple(image, path):
    if os.path.exists(path):
        os.remove(path)
    fs = FileSaver(image)
    fs.saveAsTiff(path)

def logging_broadcast(string):
	print(string)
	logging.info(string)

def check_if_files_are_present_and_equal_size(file_list):
    for path in file_list:
        if not os.path.exists(path):
            return False
        sizes = [os.path.getsize(path) / 100 for path in file_list]
        if len(set(sizes)) > 1:
            return False
    return True


###### Bigstitcher 

def get_timepoint_list_from_xml(xml_path):
    if not os.path.exists(xml_path):
        return False
    tree = ET.parse(xml_path)
    root = tree.getroot()
    timepoints = []
    for timepoint_list_elem in root.findall("./SequenceDescription/Timepoints/integerpattern"):
        if timepoint_list_elem.text == None: 
            return [0]
        timepoints = [int(x) for x in timepoint_list_elem.text.split(",")]
    return timepoints

def check_whether_timepoint_view_has_psf_in_xml(xml_path, timepoint, view):
    if not os.path.exists(xml_path):
        return False
    tree = ET.parse(xml_path)
    root = tree.getroot()
    timepoints = []
     
    if root.findall("./PointSpreadFunctions/BoundingBoxDefinition[@timepoint='%s'][@setup='%s']" % (timepoint, view)) == []:
        return False
    else:
        return True

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

def get_bounding_box_coords_from_xml(xml_path, box_name):
    if not os.path.exists(xml_path):
        return None
    tree = ET.parse(xml_path)
    root = tree.getroot()
    box_coords = None
    for box in root.findall("./BoundingBoxes/BoundingBoxDefinition[@name='%s']" % box_name):
        print(box.attrib["name"])
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

def rotate_bigstitcher_dataset(dataset_xml_path, axis_of_rotation, angle):
    logging.info("Rotating dataset around: %s for angle=%s" % (axis_of_rotation, angle))

    run_string = "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Rigid apply=[Current view transformations (appends to current transforms)] define=[Rotation around axis] same_transformation_for_all_timepoints same_transformation_for_all_angles axis_all_timepoints_channel_0_illumination_0_all_angles=%s-axis rotation_all_timepoints_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, axis_of_rotation, angle)
    logging.info("Rotating dataset run string: %s " % run_string)
    IJ.run("Apply Transformations", run_string)


def apply_transformation_bigstitcher_dataset(dataset_xml_path, affine_matrix):
    short_affine_matrix = [0 for i in range(12)]
    for i in range(len(affine_matrix) - 1):
        for j in range(len(affine_matrix[0])):
            short_affine_matrix[i * 4 + j] = affine_matrix[i][j]
    logging.info("Applying transformation with matrix %s" % short_affine_matrix)

    IJ.run("Apply Transformations", "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Affine apply=[Current view transformations (appends to current transforms)] define=Matrix same_transformation_for_all_timepoints same_transformation_for_all_angles all_timepoints_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, short_affine_matrix))

def rotate_bigstitcher_dataset_1_timepoint(dataset_xml_path, axis_of_rotation, angle, timepoint):
    logging.info("Rotating dataset around: %s for angle=%s" % (axis_of_rotation, angle))

    run_string = "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Rigid apply=[Current view transformations (appends to current transforms)] define=[Rotation around axis] same_transformation_for_all_angles axis_timepoint_%s_channel_0_illumination_0_all_angles=%s-axis rotation_timepoint_%s_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, timepoint,  axis_of_rotation, timepoint, angle)
    logging.info("Rotating dataset run string: %s " % run_string)
    IJ.run("Apply Transformations", run_string)

def apply_transformation_bigstitcher_dataset_1_timepoint(dataset_xml_path, affine_matrix, timepoint):
    short_affine_matrix = [0 for i in range(12)]
    for i in range(len(affine_matrix) - 1):
        for j in range(len(affine_matrix[0])):
            short_affine_matrix[i * 4 + j] = affine_matrix[i][j]
    logging.info("Applying transformation with matrix %s" % short_affine_matrix)
    IJ.run("Apply Transformations", "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Affine apply=[Current view transformations (appends to current transforms)] same_transformation_for_all_angles timepoint_%s_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, timepoint, short_affine_matrix))

def define_bounding_box_for_fusion(dataset_xml_path, box_coords, bounding_box_name):
    IJ.run("Define Bounding Box for Fusion", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] bounding_box=[Maximal Bounding Box spanning all transformed views] bounding_box_name=%s minimal_x=%s minimal_y=%s minimal_z=%s maximal_x=%s maximal_y=%s maximal_z=%s" % (dataset_xml_path, bounding_box_name, box_coords["x_min"], box_coords["y_min"], box_coords["z_min"], box_coords["x_max"], box_coords["y_max"], box_coords["z_max"]))

def fuse_dataset_to_tiff(dataset_xml_path, bounding_box_name, fused_xml_path, use_entropy_weighted_fusion=False):
    use_weighted = ""
    if use_entropy_weighted_fusion:
        use_weighted = " use "

    fused_dir = os.path.dirname(fused_xml_path)
    fused_tiff_names_in_dir = [os.path.basename(path) for path in get_tiffs_in_directory(fused_dir) if "fused_tp_" in path]
    already_fused_timepoints = [int(filename.split("_")[2]) for filename in fused_tiff_names_in_dir]
    all_timepoints_in_datset = get_timepoint_list_from_xml(dataset_xml_path)
    timepoints_left_to_fuse = [tp for tp in all_timepoints_in_datset if tp not in already_fused_timepoints]

    if timepoints_left_to_fuse == []:
        logging_broadcast("Found previously fused timepoints, No timepoints left to fuse!")
        return True
    if timepoints_left_to_fuse == all_timepoints_in_datset:
        timepoint_selection_string = "[All Timepoints]"
    else:
        logging_broadcast("Found previously fused timepoints, skipping them: %s" % already_fused_timepoints)
        timepoints_string = " ".join(["timepoint_%s" % tp for tp in timepoints_left_to_fuse])
        timepoint_selection_string = "[Multiple Timepoints (Select from List)] " + timepoints_string + " "
    
    fusion_string = "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=%s bounding_box=%s downsampling=1 pixel_type=[16-bit unsigned integer] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend%s produce=[Each timepoint & channel] fused_image=[Save as new XML Project (TIFF)] export_path=%s" % (dataset_xml_path, timepoint_selection_string, bounding_box_name, use_weighted, fused_xml_path)
    
    logging.info("Running fusion with following parameters\n %s" % fusion_string)
    IJ.run("Fuse dataset ...", fusion_string)
    return True
     
    

def extract_psf(dataset_xml_path, timepoint=1):
    IJ.run("Extract PSFs", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] processing_timepoint=[Timepoint %s] interest_points=beads use_corresponding remove_min_intensity_projections_from_psf psf_size_x=19 psf_size_y=19 psf_size_z=15" % (dataset_xml_path, timepoint))


def assign_psf(dataset_xml_path, from_timepoint=1):
    IJ.run("Assign PSFs", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] type=[Duplicate PSFs from other timepoint] source_timepoint=[Timepoint %s]" % (dataset_xml_path, from_timepoint)) 


def deconvolve_dataset_to_tiff(dataset_xml_path, 
                            bounding_box_name, 
                            fused_xml_path, 
                            number_deconv_iterations=4, 
                            compute_on="[CPU (Java)]", 
                            cuda_directory=None, 
                            gpu_id = None):
    #cuda_directory="/home/tema/bin/fiji-linux64/Fiji.app/lib/linux64 select_native_library_for_cudafourierconvolution=libFourierConvolutionCUDALib.so" 
    if cuda_directory is not None:
        gpu_id = "gpu_1"
        cuda_directory = "cuda_directory=" + cuda_directory + " select_native_library_for_cudafourierconvolution=libFourierConvolutionCUDALib.so "
    
    box_dims = get_bounding_box_coords_from_xml(dataset_xml_path, bounding_box_name)
    box_string = "[%s (%sx%sx%spx)]" % (bounding_box_name,
                                         abs(box_dims["x_min"] - box_dims["x_max"]) + 1,
                                         abs(box_dims["y_min"] - box_dims["y_max"]) + 1,
                                         abs(box_dims["z_min"] - box_dims["z_max"]) + 1,
                                         )

    fused_dir = os.path.dirname(fused_xml_path)
    fused_tiff_names_in_dir = [os.path.basename(path) for path in get_tiffs_in_directory(fused_dir) if "fused_tp_" in path]
    already_fused_timepoints = [int(filename.split("_")[2]) for filename in fused_tiff_names_in_dir]
    all_timepoints_in_datset = get_timepoint_list_from_xml(dataset_xml_path)
    timepoints_left_to_fuse = [tp for tp in all_timepoints_in_datset if tp not in already_fused_timepoints]

    if timepoints_left_to_fuse == []:
        logging_broadcast("Found previously fused timepoints, No timepoints left to fuse!")
        return True
    if timepoints_left_to_fuse == all_timepoints_in_datset:
        timepoint_selection_string = "[All Timepoints]"
    else:
        logging_broadcast("Found previously fused timepoints, skipping them: %s" % already_fused_timepoints)
        timepoints_string = " ".join(["timepoint_%s" % tp for tp in timepoints_left_to_fuse])
        timepoint_selection_string = "[Multiple Timepoints (Select from List)] " + timepoints_string + " "
    
    deconv_run_string = "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=%s bounding_box=%s downsampling=1 input=[Precompute Image] weight=[Precompute Image] initialize_with=[Blurred, fused image (suggested, higher compute effort)] type_of_iteration=[Efficient Bayesian - Optimization I (fast, precise)] fast_sequential_iterations osem_acceleration=1 number_of_iterations=%s use_tikhonov_regularization tikhonov_parameter=0.0060 compute=[in 512x512x512 blocks] compute_on=%s produce=[Each timepoint & channel] fused_image=[Save as new XML Project (TIFF)] %s %s export_path=%s" % (
                dataset_xml_path, 
                timepoint_selection_string,
                box_string, 
                number_deconv_iterations, 
                compute_on, 
                cuda_directory,
                gpu_id,
                fused_xml_path)
    logging.info("Running deconvolution with following parameters\n %s" % deconv_run_string)
    IJ.run("Deconvolve", deconv_run_string)


def get_transformed_PSF_from_bigstitcher_dataset(dataset_xml_path, timepoint, view):
    IJ.run("View PSFs", "select=%s process_angle=[Single angle (Select from List)] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] processing_angle=[angle %s] processing_timepoint=[Timepoint %s] display=[Averaged transformed PSF]" % (dataset_xml_path, view, timepoint))
    psf = WindowManager.getImage("Averaged transformed PSF")
    psf.hide()
    return psf


def save_transformed_psfs(dataset_xml_path, transformed_psf_dir, timepoint, num_angles):
    psf_paths = []
    for angle in range(num_angles):
        psf = get_transformed_PSF_from_bigstitcher_dataset(dataset_xml_path, timepoint, angle + 1)
        save_path = os.path.join(transformed_psf_dir, "transformed_psf_tp_%s_angle_%s.tiff" % (timepoint, angle))
        psf_paths.append(save_path)
        save_tiff_simple(psf, save_path)
    return psf_paths


def save_raw_transformed_stacks(dataset_xml_path, temp_dir_fusion, timepoint, num_angles):
    """Extract raw transformed according to registration stacks using BigStitcher Fuse method. 
    Guaranties that the stacks have been written to disk, otherwise fails and returns False.

    Args:
        dataset_xml_path (str): path to dataset xml
        temp_dir_fusion (str): path to cache directory which is supposed to be fast storage
        timepoint (int): timepoint to extract raw transformed stacks from
        num_angles (int): number of angles(views) that is going to be extracted

    Returns:
        list[str]: list of paths to raw transformed stacks
        OR
        False: if reached timeout on waiting for BigStitcher to extract the stacks
    """
    xml_path = os.path.join(temp_dir_fusion, "raw_registered_cropped.xml")
    raw_transformed_paths = [os.path.join(temp_dir_fusion, "fused_tp_%s_vs_%s.tif" % (timepoint, angle)) for angle in range(num_angles)]
    if all([os.path.exists(path) for path in raw_transformed_paths]):
        print("Found files for raw transformed stacks for the timepoint: %s Skipping raw transformed stacks creation." % timepoint) 
        return raw_transformed_paths
    
    # IJ.run("Fuse dataset ...", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] processing_timepoint=[Timepoint %s] bounding_box=embryo_cropped downsampling=1 pixel_type=[32-bit floating point] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend produce=[Each view] fused_image=[Save as new XML Project (TIFF)] export_path=%s" % (dataset_xml_path, timepoint, xml_path))
    fuse_run_string = "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] processing_timepoint=[Timepoint %s] bounding_box=embryo_cropped downsampling=1 pixel_type=[16-bit unsigned integer] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] produce=[Each view] fused_image=[Save as new XML Project (TIFF)] export_path=%s" % (dataset_xml_path, timepoint, xml_path)
    logging.info("Extracting raw transformed stack. Run string: %s" % fuse_run_string)
    IJ.run("Fuse dataset ...", fuse_run_string)


    timeout = 1200
    while check_if_files_are_present_and_equal_size(raw_transformed_paths) == False:
        time.sleep(1)
        if timeout % 10 == 0:
            logging.info("Waiting for the raw transformed stacks to be created")
        timeout -= 1
        if timeout < 0:
            logging_broadcast("ERROR: Reached timeout on waiting for raw transformed stacks to be extracted.")
            return False
    logging.info("Finished extracting raw transformed stacks.")

    return raw_transformed_paths

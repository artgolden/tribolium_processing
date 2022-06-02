# @ File(label='Fused image stack', style='file') fused_stack
# @ File(label='Raw dataset XML file', style='file') dataset_xml
# @ File(label='fused dataset XML', style='file') fused_xml
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convert

 
""" Given fused TIFF stack from 1st time point, and the original unfused dataset: Determine embryo position from max projections and fuse again cropping around the embryo. Output projections, masks and the fused cropped dataset for display."""


import shutil
from net.imagej.axis import Axes
from net.imagej.ops import Ops
from net.imagej import Dataset


from distutils.dir_util import mkpath
import math
import os
import re
import fnmatch
import json
import logging
from datetime import datetime
import traceback
import xml.etree.ElementTree as ET


from java.io import File
from ij.io import FileSaver
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin.filter import RankFilters
from ij.plugin.filter import Analyzer
from ij.plugin.filter import ParticleAnalyzer
from fiji.threshold import Auto_Threshold
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from ij.plugin import RoiRotator
from ij.plugin import RoiEnlarger
from ij.plugin import ZProjector
from ij.plugin import Slicer
from ij.plugin import StackCombiner
from ij.plugin import StackMaker
from ij.plugin import Selection
from ij.process import ByteProcessor
from ij.io import RoiDecoder
from ij.gui import PointRoi, RotatedRectRoi
from emblcmci import BleachCorrection_MH


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


def get_embryo_bounding_rectangle_and_params(max_projection, show_mask=False):

    fill_zero_pixels_with_median_intensity(max_projection)
    image_hist = max_projection.getProcessor().getHistogram()
    mean_intensity = get_mean_intensity(max_projection)
    median_intensity = get_median_intensity(max_projection)

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


def rotate_dataset(dataset_xml_path, axis_of_rotation, angle, timepoint=1):
    print("Rotating dataset around: %s for angle=%s" % (axis_of_rotation, angle))
    IJ.run("Apply Transformations", "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Rigid apply=[Current view transformations (appends to current transforms)] define=[Rotation around axis] same_transformation_for_all_angles axis_timepoint_%s_channel_0_illumination_0_all_angles=%s-axis rotation_timepoint_%s_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, timepoint,  axis_of_rotation, timepoint, angle))


def apply_transformation_dataset(dataset_xml_path, affine_matrix, timepoint=1):
    short_affine_matrix = [0 for i in range(12)]
    for i in range(len(affine_matrix) - 1):
        for j in range(len(affine_matrix[0])):
            short_affine_matrix[i * 4 + j] = affine_matrix[i][j]
    print("Applying transformation with matrix %s" % short_affine_matrix)

    IJ.run("Apply Transformations", "select=%s apply_to_angle=[All angles] apply_to_channel=[All channels] apply_to_illumination=[All illuminations] apply_to_tile=[All tiles] apply_to_timepoint=[All Timepoints] transformation=Affine apply=[Current view transformations (appends to current transforms)] same_transformation_for_all_angles timepoint_%s_channel_0_illumination_0_all_angles=%s" % (dataset_xml_path, timepoint, short_affine_matrix))


def define_bounding_box_for_fusion(dataset_xml_path, box_coords, bounding_box_name, timepoint=1):
    IJ.run("Define Bounding Box for Fusion", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] processing_timepoint=[Timepoint %s] bounding_box=[Maximal Bounding Box spanning all transformed views] bounding_box_name=%s minimal_x=%s minimal_y=%s minimal_z=%s maximal_x=%s maximal_y=%s maximal_z=%s" % (dataset_xml_path, timepoint, bounding_box_name, box_coords["x_min"], box_coords["y_min"], box_coords["z_min"], box_coords["x_max"], box_coords["y_max"], box_coords["z_max"]))


def fuse_dataset_to_display(dataset_xml_path, bounding_box_name):
    IJ.run("Fuse dataset ...", "select=%s process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] bounding_box=%s downsampling=1 pixel_type=[16-bit unsigned integer] interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend produce=[Each timepoint & channel] fused_image=[Display using ImageJ]" % (dataset_xml_path, bounding_box_name))


identity_matrix = [[0 for col in range(4)] for row in range(4)]
identity_matrix[0][0] = 1
identity_matrix[1][1] = 1
identity_matrix[2][2] = 1
identity_matrix[3][3] = 1

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


# fused_stack = IJ.openImage("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_dataset/test_rot_and_cropping/fused_tp_1_ch_0.tif")
# fused_xml_path = "/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_dataset/test_rot_and_cropping/fusion.xml"
# raw_datset_xml_path = "/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_dataset/test_rot_and_cropping/dataset__.xml"
# shutil.copy2("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_dataset/test_rot_and_cropping/dataset.xml", raw_datset_xml_path)
fused_stack = fused_stack.getAbsolutePath()
fused_xml_path = fused_xml.getAbsolutePath()
raw_dataset_xml_path = dataset_xml.getAbsolutePath()
z_projection = project_image(fused_stack, "Z", "Max")
y_projection = project_image(fused_stack, "Y", "Max")
z_bounding_roi = get_embryo_bounding_rectangle_and_params(max_projection=z_projection, show_mask=True)
y_bounding_roi = get_embryo_bounding_rectangle_and_params(max_projection=y_projection, show_mask=True)

if abs(z_bounding_roi["embryo_center_x"] - y_bounding_roi["embryo_center_x"]) > 5:
    print("Could not reliably determine embryo X center position from Z and Y projections.")

print("Embryo box params from Z projection:\n%s" % z_bounding_roi)
print("Embryo box params from Y projection:\n%s " % y_bounding_roi)

z_projection.setRoi(z_bounding_roi["bounding_roi_rect"])
y_projection.setRoi(y_bounding_roi["bounding_roi_rect"])
z_projection.show()
y_projection.show()


transformation_embryo_to_center = identity_matrix
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

apply_transformation_dataset(raw_dataset_xml_path, transformation_embryo_to_center)
embryo_crop_box = {
    "x_min" : -1 * int(z_bounding_roi["embryo_length"] / 2 + 5),
    "y_min" : -1 * int(z_bounding_roi["embryo_width"] / 2 + 5),
    "z_min" : -1 * int(y_bounding_roi["embryo_width"] / 2 + 5),
    "x_max" : int(z_bounding_roi["embryo_length"] / 2 + 5),
    "y_max" : int(z_bounding_roi["embryo_width"] / 2 + 5),
    "z_max" : int(y_bounding_roi["embryo_width"] / 2 + 5)
}

rotate_dataset(raw_dataset_xml_path, "z", round(z_bounding_roi["bounding_rect_angle"], 1))
rotate_dataset(raw_dataset_xml_path, "y", -1 * round(y_bounding_roi["bounding_rect_angle"], 1))


print(embryo_crop_box)
define_bounding_box_for_fusion(raw_dataset_xml_path, embryo_crop_box, "embryo_cropped")
fuse_dataset_to_display(raw_dataset_xml_path, "embryo_cropped")




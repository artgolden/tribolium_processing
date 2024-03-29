# @ File(label='Choose a directory with datasets', style='directory') datasets_dir
# @ File(label='Choose a file with metadata about embryo directions', style='file') metadata_file
# @ String(label='Dataset prefix', value='MGolden2022A-') dataset_name_prefix
# @ Boolean (label='Compress images?', value=true) compress_on_save
# @ Boolean (label='Use previously cropped stacks (if present)?', value=false) use_cropped_cache
# @ Boolean (label='Do histogram matching adjustment?', value=true) do_histogram_matching
# @ Integer (label='Percentage of overexposed pixels during histogram contrast adjustment', value=1) PERCENT_OVEREXPOSED_PIXELS

# Written by Artemiy Golden on Jan 2022 at AK Stelzer Group at Goethe Universitaet Frankfurt am Main
# For detailed documentation go to https://github.com/artgolden/fiji_scripts or read the README.md
# Last manual update of this line 2022.5.18 :)

from distutils.dir_util import mkpath
import math
import os
import re
import fnmatch
import json
import logging
from datetime import datetime
import traceback


from java.io import File
from ij.io import FileSaver
from ij import IJ, ImagePlus, ImageStack, WindowManager
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

#TODO:
#	What to do if embryo is off-center so much that start_plane = int(round(middle_y - 75))
# 	gives negative values? Make a better solution than just shifting the 150 plane crop.




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
			"use_manual_bounding_box": false
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
RAW_CROPPED_DIR_NAME = "(B2)-ZStacks"
CONTRAST_DIR_NAME = "(B4)-TStacks-ZN"
MONTAGE_DIR_NAME = "(B5)-TStacks-ZN-Montage"

DATASET_ERROR_FILE_NAME = "B_BRANCH_ERRORED"
DATASET_FINISHED_FILE_NAME = "B_BRANCH_FINISHED"
DATASET_ACTIVE_FILE_NAME = "B_BRANCH_ACTIVE"

MANUAL_CROP_BOX_FILE_NAME = "manual_crop_box"


DEFAULT_CROP_BOX_WIDTH = 1100
MINIMUM_CROP_BOX_WIDTH = 1000

# Percentage of overexposed pixels during histogram contrast adjustment
# PERCENT_OVEREXPOSED_PIXELS = 1


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
	logging.basicConfig(filename=os.path.join(datasets_dir, "%s-b_branch.log" % dt_string),
					 filemode='w',
					 format='%(asctime)s-%(levelname)s - %(message)s',
					 datefmt='%d-%b-%y %H:%M:%S',
					 level=logging.INFO)

	# Check metadata file for correctness
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
		if specimen_directions_in_channels == ():
			print("Error while parsing .json file: no channels found for the dataset with ID: \"%s\". Exiting." % dataset_id)
			print(EXAMPLE_JSON_METADATA_FILE)
			exit(1)
		raw_images_dir = get_raw_images_dir(datasets_dir, dataset_id)
		if not raw_images_dir:
			print("Exiting.")
			exit(1)
	logging.info("Checked metadata file. No errors found.")

	for dataset in datasets_meta["datasets"]:
		skip_the_dataset = False
		dataset_id = dataset["ID"]
		specimen_directions_in_channels = []
		for chan in range(3):
			channel = "channel_%s" % chan
			if channel in dataset:
				specimen_directions_in_channels.append(
					tuple(dataset[channel]["specimens_for_directions_1234"]))
		specimen_directions_in_channels = tuple(specimen_directions_in_channels)
		ndirections = len(specimen_directions_in_channels[0])
		
		logging.info("\n%s\nStarted processing dataset: DS%04d \n%s\n" %
					 ("#" * 100, dataset_id, "#" * 100))
		print("Started processing dataset: DS%04d" % dataset_id)
		raw_images_dir = get_raw_images_dir(datasets_dir, dataset_id)
		root_dataset_dir = os.path.split(raw_images_dir)[0]
		if os.path.exists(os.path.join(root_dataset_dir, DATASET_FINISHED_FILE_NAME)):
			logging.info(
				"Found %s file. Dataset DS%04d already processed, skipping." % (DATASET_FINISHED_FILE_NAME, dataset_id))
			print("Found %s file. Dataset DS%04d already processed, skipping." %
				  (DATASET_FINISHED_FILE_NAME, dataset_id))
			continue
		if os.path.exists(os.path.join(root_dataset_dir, DATASET_ERROR_FILE_NAME)):
			logging.info(
				"Found %s file. Dataset DS%04d errored while previous processing, skipping." % (DATASET_ERROR_FILE_NAME, dataset_id))
			print("Found %s file. Dataset DS%04d errored while previous processing, skipping." %
				  (DATASET_ERROR_FILE_NAME, dataset_id))
			continue
		if os.path.exists(os.path.join(root_dataset_dir, DATASET_ACTIVE_FILE_NAME)):
			logging.info(
				"Found %s file. Perhaps the dataset DS%04d is currently being processed by other Fiji instance, skipping." % (DATASET_ACTIVE_FILE_NAME, dataset_id))
			print("Found %s file. Perhaps the dataset DS%04d is currently being processed by other Fiji instance, skipping." %
				  (DATASET_ACTIVE_FILE_NAME, dataset_id))
			continue
		open(os.path.join(root_dataset_dir, DATASET_ACTIVE_FILE_NAME), 'a').close()

		logging.info("\tArranging raw image files")
		try:
			move_files(
				raw_images_dir, specimen_directions_in_channels, dataset_id, dataset_name_prefix)
		except ValueError as e:
			logging.info(
				"Error while moving files for the dataset:\"%s\", skipping the dataset. Error:\n %s" % (dataset_id, e))
			print("Error while moving files for the dataset:\"%s\", skipping the dataset." % dataset_id)
			print(e)
			continue
		logging.info("\tFiles arranged, starting processing.")

		# We have to make sure that all directions in all channels have the same crop box dimensions.
		dataset_minimal_crop_box_dims = (0, 0)
		channel = 1
		chan_dir_name = "CH%04d" % channel
		raw_images_direction_dirs = make_directions_dirs(
			os.path.join(raw_images_dir, chan_dir_name), ndirections)
		meta_dir = os.path.join(root_dataset_dir, METADATA_DIR_NAME)
		meta_d_dirs = make_directions_dirs(os.path.join(meta_dir, chan_dir_name), ndirections)
		tstack_dataset_dirs = make_directions_dirs(
					os.path.join(root_dataset_dir, TSTACKS_DIR_NAME, chan_dir_name), ndirections)
		raw_cropped_dirs = make_directions_dirs(os.path.join(
			root_dataset_dir, RAW_CROPPED_DIR_NAME, chan_dir_name), ndirections)
		# Loop one time over all dimensions in channel 1 to determine the dataset_maximal_crop_box_width
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
			mproj_stack_file_name = get_tiff_name_from_dir(raw_dir)
			mproj_stack_file_name.plane = "(ZM)"
			mproj_stack_file_name.time_point = "(TS)"
			mproj_stack_path = os.path.join(
				tstack_backup_dir, mproj_stack_file_name.get_name())
			if not os.path.exists(mproj_stack_path):
				max_proj_stack = make_max_Z_projections_for_folder(raw_dir)
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
				max_time_proj = project_a_stack(max_proj_stack)
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
					max_time_proj, m_dir, dataset, dataset_minimal_crop_box_dims, use_dataset_box_dims=False)
			except Exception as e:
				logging.info(
					"ERROR: Encountered an exception while trying to create a crop template. Skipping the dataset. Exception:\n%s" % e)
				print("ERROR: Encountered an exception while trying to create a crop template. Skipping the dataset. Exception:\n%s" % e)
				skip_the_dataset = True
				open(os.path.join(root_dataset_dir, DATASET_ERROR_FILE_NAME), 'a').close()
				os.remove(os.path.join(root_dataset_dir, DATASET_ACTIVE_FILE_NAME))
				break
			fs = FileSaver(cropped_max_time_proj)
			fs.saveAsTiff(os.path.join(
				m_dir, "cropped_max_time_proj.tif"))

		# Main loop over the channels and directions for the dataset
		for channel, _ in enumerate(specimen_directions_in_channels, start=1):
			if skip_the_dataset == True:
				break
			chan_dir_name = "CH%04d" % channel
			raw_images_direction_dirs = make_directions_dirs(
				os.path.join(raw_images_dir, chan_dir_name), ndirections)
			root_dataset_dir = os.path.split(raw_images_dir)[0]
			meta_dir = os.path.join(root_dataset_dir, METADATA_DIR_NAME)
			meta_d_dirs = make_directions_dirs(os.path.join(meta_dir, chan_dir_name), ndirections)
			tstack_dataset_dirs = make_directions_dirs(
				os.path.join(root_dataset_dir, TSTACKS_DIR_NAME, chan_dir_name), ndirections)
			raw_cropped_dirs = make_directions_dirs(os.path.join(
				root_dataset_dir, RAW_CROPPED_DIR_NAME, chan_dir_name), ndirections)
			contrast_dirs = make_directions_dirs(os.path.join(
				root_dataset_dir, CONTRAST_DIR_NAME, chan_dir_name), ndirections)

			for direction, raw_dir, tstack_dir, m_dir, raw_cropped_dir in zip(range(1, ndirections + 1), raw_images_direction_dirs, tstack_dataset_dirs, meta_d_dirs, raw_cropped_dirs):
				if skip_the_dataset == True:
					break
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
				mproj_stack_file_name = get_tiff_name_from_dir(raw_dir)
				mproj_stack_file_name.plane = "(ZM)"
				mproj_stack_file_name.time_point = "(TS)"
				mproj_stack_path = os.path.join(
					tstack_backup_dir, mproj_stack_file_name.get_name())
				if not os.path.exists(mproj_stack_path):
					max_proj_stack = make_max_Z_projections_for_folder(raw_dir)
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
					max_time_proj = project_a_stack(max_proj_stack)
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
						max_time_proj, m_dir, dataset, dataset_minimal_crop_box_dims, use_dataset_box_dims=True)
				except Exception as e:
					traceback.print_exc()
					logging.info(
						"ERROR: Encountered an exception while trying to create a crop template. Skipping the dataset. Exception:\n%s" % e)
					print("ERROR: Encountered an exception while trying to create a crop template. Skipping the dataset. Exception:\n%s" % e)
					skip_the_dataset = True
					open(os.path.join(root_dataset_dir, DATASET_ERROR_FILE_NAME), 'a').close()
					os.remove(os.path.join(root_dataset_dir, DATASET_ACTIVE_FILE_NAME))
					break
				fs = FileSaver(cropped_max_time_proj)
				fs.saveAsTiff(os.path.join(
					m_dir, "cropped_max_time_proj.tif"))

				# Cropping max projections

				cropped_tstack_file_name = get_tiff_name_from_dir(
					raw_dir)
				cropped_tstack_file_name.time_point = "(TS)"
				cropped_tstack_file_name.plane = "(ZM)"
				if not use_cropped_cache or not os.path.exists(os.path.join(tstack_dir, cropped_tstack_file_name.get_name())):
					logging.info("\tChannel: %s Direction: %s Cropping a stack of max projections." % (
						channel, direction))
					cropped_tstack = crop_stack_by_template(
						max_proj_stack, crop_template, dataset)
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
				if use_cropped_cache and all_stacks_done:
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
							raw_stack, crop_template, dataset)
						if direction == 0:
							logging.info("\tChannel: %s Direction: %s Finding which planes to keep in raw image stacks." % (
								channel, direction))
							if "planes_to_keep_per_direction" in dataset:
								manual_planes_to_keep = (dataset["planes_to_keep_per_direction"][direction - 1]
														["start"], dataset["planes_to_keep_per_direction"][direction - 1]["end"])
							else:
								manual_planes_to_keep = None
							stack_for_finding_planes = IJ.openImage(raw_stack_files[0])
							IJ.run(stack_for_finding_planes, "Properties...",
												"frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=4.0000")
							stack_for_finding_planes = crop_stack_by_template(
															stack_for_finding_planes, crop_template, dataset)
							try:
								planes_kept = find_planes_to_keep(stack_for_finding_planes, m_dir, manual_planes_to_keep)
							except Exception as e:
								logging.info("\tChannel: %s Direction: %s Encountered an exception while trying to find which planes to keep. Skipping the dataset. Exception: \n%s" % (
									channel, direction, e))
								print("\tChannel: %s Direction: %s Encountered an exception while trying to find which planes to keep. Skipping the dataset. Exception: \n%s" % (
									channel, direction, e))
								skip_the_dataset = True
								open(os.path.join(root_dataset_dir, DATASET_ERROR_FILE_NAME), 'a').close()
								os.remove(os.path.join(root_dataset_dir, DATASET_ACTIVE_FILE_NAME))
								print("Encountered an exception while trying to find which planes to keep. Skipping the dataset.")
								break
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
					tstack_dir, get_tiff_name_from_dir(tstack_dir).get_name())
				imp_stack = IJ.openImage(stack_path)
				if direction == 0:
					montage_stack = imp_stack.getStack()
					montage_stack_name = get_tiff_name_from_dir(tstack_dir)
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
			if do_histogram_matching:
				match_histograms_stack(adj_montage_stack)
				adj_montage_stack_name.plane = "(ZH)"
			adj_montage_stack, histogram_thresholds = threshold_histogram_stack(adj_montage_stack)
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
				adj_stack_name = get_tiff_name_from_dir(tstack_dir)
				adj_stack_name.plane = "(ZA)"
				if do_histogram_matching:
					adj_stack_name.plane = "(ZH)"
			
				save_tiff(adj_max_proj_stack, os.path.join(contr_dir, adj_stack_name.get_name()))
			
			# Save adjusted vertical montages

			logging.info("\tChannel: %s Creating adjusted vertical montages from max projections." % (
                            channel))
			vert_montage_stack = ImagePlus()
			vert_montage_stack_name = None
			for direction, tstack_dir in enumerate(tstack_dataset_dirs):
				stack_path = os.path.join(
					tstack_dir, get_tiff_name_from_dir(tstack_dir).get_name())
				imp_stack = IJ.openImage(stack_path)
				if direction == 0:
					vert_montage_stack = imp_stack.getStack()
					vert_montage_stack_name = get_tiff_name_from_dir(tstack_dir)
				if direction > 0:
					vert_montage_stack = StackCombiner.combineVertically(StackCombiner(), vert_montage_stack, imp_stack.getStack())
			vert_montage_stack_name.direction = "(AY)"
			vert_montage_stack_name.plane = "(ZA)"
			vert_montage_stack = ImagePlus("montage", vert_montage_stack)
			if do_histogram_matching:
				match_histograms_stack(vert_montage_stack)
				vert_montage_stack_name.plane = "(ZH)"
			vert_montage_stack, _ = threshold_histogram_stack(vert_montage_stack)
			montage_dir = os.path.join(root_dataset_dir, MONTAGE_DIR_NAME)
			if not os.path.exists(montage_dir):
				mkpath(montage_dir)
			save_tiff(vert_montage_stack, os.path.join(montage_dir, vert_montage_stack_name.get_name()))

		if skip_the_dataset == True:
			logging.info("Had to skip the dataset DS%04d." % dataset_id)
			continue
		open(os.path.join(root_dataset_dir, DATASET_FINISHED_FILE_NAME), 'a').close()
		os.remove(os.path.join(root_dataset_dir, DATASET_ACTIVE_FILE_NAME))
		logging.info("Finished processing dataset DS%04d successfully." % dataset_id)
		print("Finished processing dataset DS%04d successfully." % dataset_id)

	logging.info("Finished processing all datasets.")
	print("Finished processing all datasets.")


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


def make_directions_dirs(input_dir, ndirections):
	direction_dirs = []
	for i in range(1, ndirections + 1):
			new_dir = os.path.join(input_dir, "DR" + str(i).zfill(4))
			direction_dirs.append(new_dir)
			if not os.path.exists(new_dir):
				mkpath(new_dir)

	return direction_dirs


def get_tiff_name_from_dir(input_dir):
	file_name = None
	for fname in os.listdir(input_dir):
		if fname.lower().endswith(".tif"):
			file_name = fname
			break
	if file_name == None:
		raise Exception("Did not find any TIFF files in a directory: %s" % input_dir)

	file_name = FredericFile(file_name)
	return file_name


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
		direction_dirs = make_directions_dirs(chan_dir, ndirections)
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


def is_specimen_input_valid(specimens_per_direction):
	if not isinstance(specimens_per_direction, tuple):
		return False
	for i in specimens_per_direction:
		if not isinstance(i, int):
			return False
		if i not in range(1000):
			return False
	return True


def is_dataset_ID_input_valid(dataset_id):
	if not isinstance(dataset_id, int):
		return False
	if dataset_id in range(10000):
		return True
	return False


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


def sort_tiff_list_by_timepoint(tiff_list):
	pass


def project_a_stack(stack):
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


def make_max_Z_projections_for_folder(input_dir):
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
		max_proj = project_a_stack(stack)
		stack_list.append(max_proj.getProcessor())
	max_proj_stack = ImageStack(img_for_dims.width, img_for_dims.height)
	for slice in stack_list:
		max_proj_stack.addSlice(None, slice)
	max_proj_stack = ImagePlus("max_proj", max_proj_stack)
	max_proj_stack = reset_img_properties(max_proj_stack, voxel_depth=1)
	return max_proj_stack


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


def create_crop_template(max_time_projection, meta_dir, dataset, dataset_minimal_crop_box_dims, use_dataset_box_dims):
	"""Crop max_time_projection, create crop template .roi object, save it in meta_dir.

	Args:
		max_time_projection (ImagePlus): a single image of max projection of a time stack of max projections
		meta_dir (str): full path to metadata dir
		dataset (dict): metadata about the dataset to extract information about embryo orientation
		dataset_maximal_crop_box_dims (int, int): width, heigth
		use_dataset_box_dims (bool):

	Returns:
		Roi: bounding box around the embryo with proper rotation and size
		ImagePlus: cropped max time projection
		(int, int): updated dataset maximal crop box dims
	"""
	updated_dataset_maximal_crop_box_dims = dataset_minimal_crop_box_dims
	imp = max_time_projection

	if dataset["use_manual_bounding_box"] == True:
		logging.info("\tUsing manualy specified crop box.")
		manual_roi = RoiDecoder(os.path.join(meta_dir, "%s.roi" % MANUAL_CROP_BOX_FILE_NAME)).getRoi()
		rot_angle = round(get_polygon_roi_angle(manual_roi), ndigits=1)
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
		bt = mask.getProcessor().convertToByteProcessor()
		for i in range(40):
			bt.erode(2, 0)
		for i in range(40):
			bt.dilate(2, 0)
		for i in range(40):
			bt.erode(3, 0)
		for i in range(40):
			bt.dilate(3, 0)
		mask.setProcessor(bt)
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
	if dataset["use_manual_bounding_box"] == False:
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
	if dataset["head_direction"] == "left":
		IJ.run(cropped_max_time_proj, "Rotate 90 Degrees Right", "")

		# So that illumination is comming from the right.
		# Right now everything is based on that illumination on the images from mDSLM is always coming from the top.
		IJ.run(cropped_max_time_proj, "Flip Horizontally", "")
	else:
		IJ.run(cropped_max_time_proj, "Rotate 90 Degrees Left", "")
	

	return (crop_template, cropped_max_time_proj, updated_dataset_maximal_crop_box_dims)


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


def midpoint(x1, y1, x2, y2):
	x1, y1, x2, y2 = float(x1), float(y1), float(x2), float(y2)
	return ((x1 + x2) / 2, (y1 + y2) / 2)


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


def crop_stack_by_template(stack, crop_template, dataset):
	"""Crop stack by provided PolygonRoi which should be a rotated rectangle around the embryo.
	It also flips the ebryo accroding to metadata about embryo facing direction, so that illumination is from the left.

	Args:
		stack (ImagePlus): a stack of image planes
		crop_template (PolygonRoi): rotated bounding rectangle around the embryo
		dataset (dict): metadata about the dataset to extract information about embryo orientation

	Returns:
		ImagePlus: cropped stack
	"""

	stack.setRoi(crop_template)
	cropped_stack = stack.crop("stack")
	stack = None
	# cropped_stack.show()
	IJ.run(cropped_stack, "Select All", "")

	IJ.run(cropped_stack, "Rotate... ",
			"angle=%s grid=1 interpolation=Bilinear stack" % round(get_polygon_roi_angle(crop_template), ndigits=1))
	final_center_x = cropped_stack.getWidth() / 2
	final_center_y = cropped_stack.getHeight() / 2
	box_width, box_height = get_rotated_rect_polygon_roi_dims(crop_template)
	IJ.run(cropped_stack, "Specify...", "width=%s height=%s x=%s y=%s centered" % (box_width, box_height, final_center_x, final_center_y))
	cropped_stack_resized = cropped_stack.crop("stack")

	if dataset["head_direction"] == "left":
		IJ.run(cropped_stack_resized, "Rotate 90 Degrees Right", "")
		IJ.run(cropped_stack_resized, "Flip Horizontally", "stack")
	else:
		IJ.run(cropped_stack_resized, "Rotate 90 Degrees Left", "")
	IJ.run(cropped_stack_resized, "Select None", "")
	return cropped_stack_resized


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

		max_proj_reslice = project_a_stack(resliced)
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

	for_user_assessment = project_a_stack(for_user_assessment)
	fs = FileSaver(for_user_assessment)
	fs.saveAsTiff(os.path.join(
		meta_dir, "Y_projected_raw_stack_for_assessment_of_plane_selection.tif"))

	return (start_plane, end_plane)


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


def get_histogram_thresholds(image_processor):
	"""Get upper and lower thresholds for the histogram of the image for contrast adjustment. 
	Lower histogram is determined by Triangle thresholding algorithm.
	Upper histogram is determined in such a way that globally determined persentage of pixels in the image after applying lower threshold will be overexposed.

	Args:
		image_processor (ImageProcessor): image to determine histogram thresholds from

	Returns:
		(int, int): a tuple of lower and upper thresholds
	"""
	hist = image_processor.getHistogram()
	lower_threshold = Auto_Threshold.Triangle(hist)
	num_overexposed_pixels = 0
	upper_threshold = 65535
	sum_elem_above_threshold = sum(hist[i] for i in range(lower_threshold, len(hist)))

	# making it so 1% of all pixels will be overexposed
	num_overexposed_pixels_threshold = float(PERCENT_OVEREXPOSED_PIXELS) / 100 * sum_elem_above_threshold
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


def threshold_histogram_stack(imp_stack):
	"""Adjust contrast of the image stack by removing pixels based on the lower and upper intensity thresholds. Intensity thresholds are determined from the middle plane in the stack by get_histogram_thresholds function. 

	Args:
		imp_stack (ImagePlus): image stack to adjust

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
	histogram_thresholds = get_histogram_thresholds(stack.getProcessor(middle_plane))
	logging.info("Adjusted stack with thresholds: %s, %s." % (
            histogram_thresholds[0], histogram_thresholds[1]))
	for i in range(1, nplanes + 1):
		imp2 = ImagePlus("slice_iterator", stack.getProcessor(i))
		auto_contrast_by_histogram_thresholds(imp2, histogram_thresholds)
		adjusted_stack.append(imp2)
	adjusted_stack = ImagePlus(
		"Adjusted contrast", ImageStack.create(adjusted_stack))
	return (adjusted_stack, histogram_thresholds)


def save_tiff(image, path):
	if os.path.exists(path):
		os.remove(path)
	if compress_on_save == False:
		fs = FileSaver(image)
		fs.saveAsTiff(path)
	else:
		IJ.run(image, "Bio-Formats Exporter", "save=%s export compression=zlib" % path)


def match_histograms_stack(stack):
	"""Adjust histrograms for the provided image stack. Mutates provided stack. Resulting histograms for each plane have the same area
	under the curve and look equally bright.

	Args:
		stack (ImagePlus): a stack of images
	"""
	bleach_obj = BleachCorrection_MH(stack)
	bleach_obj.doCorrection()


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

if __name__ in ['__builtin__', '__main__']:
	process_datasets(datasets_dir, metadata_file, dataset_name_prefix)

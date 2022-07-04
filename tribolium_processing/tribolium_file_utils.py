from copy import deepcopy
import os
import re


class FredericFile:
	"""
	File naming for light-sheet image files. by Frederic Strobl 
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
	additional_info = ""

	def __init__(self, file_name):
		split_by = "-DS|TP|DR|CH|PL|\."
		name_parts = re.split(split_by, file_name)
		
		if len(name_parts) == 7:
			self.dataset_name, self.dataset_id, self.time_point, self.direction, self.channel, plane_and_info, self.extension = name_parts
			plane_and_info_list = plane_and_info.split("_", 1) # Split by first occurance only
			if len(plane_and_info_list) == 1:
				self.plane = plane_and_info_list[0]
			else:
				self.plane, self.additional_info = plane_and_info_list
		else:
			raise Exception(
				"Image file name is improperly formatted! Check documentation inside the script. Expected 7 parts after splitting by %s" % split_by)

		self.extension.lower()

	def get_name(self):
		if self.additional_info != "":
			self.additional_info = "_" + self.additional_info
		return "%s-DS%sTP%sDR%sCH%sPL%s%s.%s" % (self.dataset_name,
										 self.dataset_id,
										 self.time_point,
										 self.direction,
										 self.channel,
										 self.plane,
										 self.additional_info,
										 self.extension)
	def get_name_without_extension(self):
		if self.additional_info != "":
			self.additional_info = "_" + self.additional_info
		return "%s-DS%sTP%sDR%sCH%sPL%s%s" % (self.dataset_name,
										 self.dataset_id,
										 self.time_point,
										 self.direction,
										 self.channel,
										 self.plane,
										 self.additional_info)

	def set_param(self, dataset_name=None, dataset_id=None, time_point=None, direction=None, channel=None, plane=None, extension=None, additional_info=None):
		if dataset_name is not None:
			self.dataset_name = dataset_name

		if dataset_id is not None:
			if isinstance(dataset_id, int):
				dataset_id = "%04d" % dataset_id
			self.dataset_id = dataset_id

		if time_point is not None:
			if isinstance(time_point, int):
				time_point = "%04d" % time_point
			self.time_point = time_point

		if direction is not None:
			if isinstance(direction, int):
				direction = "%04d" % direction
			self.direction = direction

		if channel is not None:
			if isinstance(channel, int):
				channel = "%04d" % channel
			self.channel = channel

		if plane is not None:
			if isinstance(plane, int):
				plane = "%04d" % plane
			self.plane = plane

		if extension is not None:
			self.extension = extension
			
		if additional_info is not None:
			self.additional_info = additional_info

	def get_modified_name(self, dataset_name=None, dataset_id=None, time_point=None, direction=None, channel=None, plane=None, extension=None, additional_info=None):
		changed = deepcopy(self)
		changed.set_param(dataset_name, dataset_id, time_point, direction, channel, plane, extension, additional_info)
		return changed.get_name()

	def get_modified_name_without_extension(self, dataset_name=None, dataset_id=None, time_point=None, direction=None, channel=None, plane=None, extension=None, additional_info=None):
		changed = deepcopy(self)
		changed.set_param(dataset_name, dataset_id, time_point, direction, channel, plane, extension, additional_info)
		return changed.get_name_without_extension()



###### File managing

def get_tiffs_in_directory(directory):
	"""Get all TIFF file paths in a directory. Subdirectories are not searched.

	Args:
		directory (str): full path to directory

	Returns:
		str[]: list of full paths to tiff files
	"""
	file_names = []
	for fname in os.listdir(directory):
		if fname.lower().endswith(".tif") or fname.lower().endswith(".tiff"):
			file_names.append(os.path.join(directory, fname))
	file_names = sorted(file_names)
	return file_names

def get_any_tiff_name_from_dir(input_dir):
	file_name = None
	for fname in os.listdir(input_dir):
		if fname.lower().endswith(".tif") or fname.lower().endswith(".tiff"):
			file_name = fname
			break
	if file_name == None:
		raise Exception("Did not find any TIFF files in a directory: %s" % input_dir)

	file_name = FredericFile(file_name)
	return file_name
	
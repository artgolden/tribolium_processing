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
	def get_name_without_extension(self):
		return "%s-DS%sTP%sDR%sCH%sPL%s" % (self.dataset_name,
										 self.dataset_id,
										 self.time_point,
										 self.direction,
										 self.channel,
										 self.plane)



###### File managing

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

def get_any_tiff_name_from_dir(input_dir):
	file_name = None
	for fname in os.listdir(input_dir):
		if fname.lower().endswith(".tif"):
			file_name = fname
			break
	if file_name == None:
		raise Exception("Did not find any TIFF files in a directory: %s" % input_dir)

	file_name = FredericFile(file_name)
	return file_name
	
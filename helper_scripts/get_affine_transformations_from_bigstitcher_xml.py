# File(label='BigStitcher dataset XML file', style='file') dataset_xml

import xml.etree.ElementTree as ET

from ij import IJ

def get_tranformations_from_xml_file(xml_path, timepoint=23, view=1 ):
	transformations = []
	tree = ET.parse(xml_path)
	root = tree.getroot()
	for view_transform in reversed(root.findall("./ViewRegistrations/ViewRegistration[@timepoint='%s'][@setup='%s']/ViewTransform/affine" % (timepoint, view))):
		transformation = [float(i) for i in view_transform.text.split(" ")]
		transformations.append(transformation)
	return transformations


print(get_tranformations_from_xml_file("/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/(P0)-ZStacks-Raw/CH0001/bigstitcher_dataset/MGolden2022A-DS0016TPall_DRall_CH0001PL(ZS).xml"))
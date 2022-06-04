# @ File(label='dataset XML file', style='file') xml_file
# @ String(label='Box name', value='max_box') box_name
# @ Integer (label='dummy variable', value=1) dummy

"""Get dictionary with coordinates from the BigStitcher dataset XML file"""

import xml.etree.ElementTree as ET

def get_bounding_box_coords_from_xml(xml_path, box_name):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    for box in root.iter('BoundingBoxes'):
         if box[0].attrib["name"] == box_name:
             print(box[0].attrib["name"])
             min = box[0].find("min").text.split(" ")
             max = box[0].find("max").text.split(" ")
             box_coords = {
                 "x_min" : int(min[0]),
                 "y_min" : int(min[1]),
                 "z_min" : int(min[2]),
                 "x_max" : int(max[0]),
                 "y_max" : int(max[1]),
                 "z_max" : int(max[2])
             }
    return box_coords

xml_file = xml_file.getAbsolutePath()
# xml_file = "/run/user/1000/gvfs/smb-share:server=bugcube.physbio.uni-frankfurt.de,share=bugcube/TemporaryData/Artemiy/test_fusion_thresholding/MGolden2022A_DS0039_Sqh5_old_DS0019_manual_crop/dataset.xml"

print get_bounding_box_coords_from_xml(xml_path=xml_file, box_name=box_name)


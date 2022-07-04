# @ File(label='dataset XML file', style='file') xml_file
# @ Integer (label='dummy variable', value=1) dummy

"""Get dictionary with coordinates from the BigStitcher dataset XML file"""

import os
import xml.etree.ElementTree as ET



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



xml_file = xml_file.getAbsolutePath()

print(check_whether_timepoint_view_has_psf_in_xml(xml_file, 48, 3))


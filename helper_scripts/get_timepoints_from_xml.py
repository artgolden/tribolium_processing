# @ File(label='dataset XML file', style='file') xml_file
# @ Integer (label='dummy variable', value=1) dummy

"""Get dictionary with coordinates from the BigStitcher dataset XML file"""

import os
import xml.etree.ElementTree as ET

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

xml_file = xml_file.getAbsolutePath()

print(get_timepoint_list_from_xml(xml_path=xml_file))


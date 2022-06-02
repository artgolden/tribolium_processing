# @ File(label='dataset XML file', style='file') xml_file
# @ Integer (label='dummy variable', value=1) dummy

"""Get matrix of transformation from fusion coordinates back to raw registered views from the BigStitcher dataset XML file"""

import xml.etree.ElementTree as ET

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

xml_file = xml_file.getAbsolutePath()

print(get_fusion_tranformation_from_xml_file(xml_path=xml_file))


import xml.etree.ElementTree as ET
import glob

# Define the xlink namespace
namespaces = {'xlink': 'http://www.w3.org/1999/xlink'}

# Open the output file
with open('output.txt', 'w') as f:
    # Loop over all SVG files in the current directory
    for filename in glob.glob('*.svg'):
        # Parse the SVG file
        tree = ET.parse(filename)
        root = tree.getroot()

        # Find all elements with an xlink:href attribute
        elements = root.findall('.//*[@xlink:href]', namespaces=namespaces)

        # Extract the xlink:href attributes
        hrefs = [element.get('{http://www.w3.org/1999/xlink}href') for element in elements]

        # Write the filename and the xlink:href attributes to the output file
        f.write(f'File: {filename}\n')
        for href in hrefs:
            f.write(f'{href}\n')
        f.write('\n')
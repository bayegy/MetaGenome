#!/usr/bin/env python3.8
'''
20201102,tangmaomao wrote this script.
You need to copy this script file into the "result" directory, and write the specified arguments, and run the script in that directory.
'''
from MetaGenome.visualizeAll import VisualizeAll
# the "Category1" may be modified
v = VisualizeAll("mapping_file.txt", "Category1", "./")
# if the result was form "based on assembly", you need to write the argument base_on_assembly=True
v.visualize(base_on_assembly=False)

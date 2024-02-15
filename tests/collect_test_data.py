"""
Script to setup test data for specified dates. This script will replace the specified test data if it exists already, after a prompt
"""

import os
import shutil

from mhdpy.analysis.standard_import import *

CHECK_FOR_EXISTING = False

def copy_files(src_folder, dest_folder):
    for root, dirs, files in os.walk(src_folder):
        subdir = os.path.relpath(root, src_folder)

        for file in files:
            fp_old = os.path.join(root, file)
            fp_new = os.path.join(dest_folder, subdir, file)

            if not os.path.exists(os.path.dirname(fp_new)):
                os.makedirs(os.path.dirname(fp_new))
            if os.path.exists(fp_new):
                if CHECK_FOR_EXISTING:
                    print('File {} already exists in test_data_folder. Skipping...'.format(fp_new))
                    continue
                else:
                    print('File {} already exists in test_data_folder. Overwriting...'.format(fp_new))
                    os.remove(fp_new)
            shutil.copy(fp_old, fp_new)
            print('Copied {} to {}'.format(fp_old, fp_new))

test_data_folder = 'test_data'
if not os.path.exists(test_data_folder):
    os.mkdir(test_data_folder)

test_data_folder_modeling = 'test_data_modeling'
if not os.path.exists(test_data_folder_modeling):
    os.mkdir(test_data_folder_modeling)

test_data_folder_final = 'test_data_final'
if not os.path.exists(test_data_folder_final):
    os.mkdir(test_data_folder_final)

input_data_folder = os.path.join(REPO_DIR, 'experiment', 'data', 'proc_data')
copy_files(input_data_folder, test_data_folder)

modeling_data_folder = os.path.join(REPO_DIR, 'modeling', 'dataset', 'output')
copy_files(modeling_data_folder, test_data_folder_modeling)

final_data_folder = os.path.join(REPO_DIR, 'experiment', 'figure_panels', 'output')
copy_files(final_data_folder, test_data_folder_final)
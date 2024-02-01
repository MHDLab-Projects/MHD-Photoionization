"""
Script to setup test data for specified dates. This script will replace the specified test data if it exists already, after a prompt
"""

import os
import shutil

from mhdpy.analysis.standard_import import *

CHECK_FOR_EXISTING = False

test_data_folder = 'test_data'
if not os.path.exists(test_data_folder): os.mkdir(test_data_folder)

input_data_folder = pjoin(REPO_DIR, 'experiment','data','proc_data')

# copy all files in input_data_folder to test_data_folder recursively, overwriting if necessary

for root, dirs, files in os.walk(input_data_folder):
    subdir = os.path.relpath(root, input_data_folder)

    if subdir == 'lecroy':
        # skip lecroy folder for now
        continue

    for file in files:
        fp_old = pjoin(root, file)

        fp_new = pjoin(test_data_folder, subdir, file)

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

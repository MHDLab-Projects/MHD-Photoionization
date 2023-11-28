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

shutil.copy(pjoin(input_data_folder, 'ds_absem.cdf'), pjoin(test_data_folder, 'ds_absem.cdf'))

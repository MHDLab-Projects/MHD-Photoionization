"""
Script to extract line and beam profiles from CFD simulation data with new Yeq fields
TODO: replace old extract_line_profiles.py with this script. Waiting on K=1% data. 
"""

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os
import pint_pandas
import pandas as pd

from mhdpy.analysis.standard_import import *
from pi_paper_utils.fileio import cfd_fp_dict, cfd_all_fields
# from mhdpy.pyvista_utils import AxiMesh, pv_to_unstack_xr, downsel_arrays, pv_to_xr
from mhdpy.pyvista_utils.axi import AxiInterpolator

create_standard_folders()

#%%[markdown]

# Demo

#%%

fp = cfd_fp_dict['0.8_1.0']

mesh = pv.read(fp)

mesh.set_active_scalars('CO2')

#%%

spacing = 0.001

bounds = (0.2, 0.6, -0.05, 0.05, 0.0, 0.0)

l_x = bounds[1] - bounds[0]
l_y = bounds[3] - bounds[2]

n_x = int(l_x/spacing)
n_y = int(l_y/spacing)

grid = pv.ImageData()
grid.origin = (bounds[0], bounds[2], 0)
grid.spacing = (spacing, spacing, 0)
grid.dimensions = (n_x, n_y, 1)

#%%


#%%

ai = AxiInterpolator(mesh, var_names=cfd_all_fields)

f_out = ai(grid.points)

#%%

df = pd.DataFrame(f_out, columns=ai.names)

df['pos_x'] = grid.points[:,0]
df['pos_y'] = grid.points[:,1]

df = df.set_index(['pos_x', 'pos_y'])

df.to_csv('output/mdot0130_phi080_K100.csv')

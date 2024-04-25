"""
Script to extract line and beam profiles from CFD simulation data with new Yeq fields
TODO: replace old extract_line_profiles.py with this script. Waiting on K=1% data. 
"""


#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from mhdpy.fileio import gen_path

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pint_pandas
from pi_paper_utils.fileio import cfd_fp_dict


soi = ['K', 'Kp', 'em', 'OH', 'OHm', 'KOH', 'O2', 'H2O', 'N2', 'CO2']
# soi_Yeq = ['Yeq_K', 'Yeq_K+', 'Yeq_e-', 'Yeq_OH', 'Yeq_OH-', 'Yeq_KOH', 'Yeq_K2CO3']
# additional = ['T', 'p', 'rho']
# all_fields = [*soi, *soi_Yeq, *additional]



# soi = ['K', 'K_p', 'e_m', 'OH', 'OH_m', 'KOH', 'O2', 'H2O', 'N2', 'CO2']
soi_Yeq = ['Yeq_K', 'Yeq_K+', 'Yeq_e-', 'Yeq_OH', 'Yeq_OH-', 'Yeq_KOH', 'Yeq_K2CO3', 'Yeq_KO']
additional = ['T', 'p', 'rho']
all_fields = [*soi, *soi_Yeq, *additional]



# %%

from mhdpy.pyvista_utils import AxiMesh, pv_to_unstack_xr, downsel_arrays, pv_to_xr



#%%[markdown]

# Demo

#%%

# fp = pjoin(results_dir, 'medium', 'mdot0130_phi080_K100', 'frontCyl_chem1.vtk')
# fp = pjoin(results_dir, 'medium', 'mdot0130_phi080_K010', 'frontCyl_chem1.vtk')
fp = cfd_fp_dict['0.8_1.0']

mesh = pv.read(fp)

mesh.set_active_scalars('CO2')

#%%

spacing = 0.001
# bounds = mesh.GetBounds()

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

grid.points

#%%

#TODO: add to mhdpy
sys.path.append(pjoin(REPO_DIR, 'modeling', 'cfd'))
from pv_axi_utils import AxiInterpolator

#%%

ai = AxiInterpolator(mesh, var_names = all_fields)

# ai returns a 2D array of the fields at the points in the grid.
f_out = ai(grid.points)

#%%

# Create a 2D array where the first dimension corresponds to the 'pos_x' and 'pos_y' coordinates,
# and the second dimension corresponds to the 'T' and 'K' variables.
# data = np.stack([f_out[:,i] for i in range(f_out.shape[1])], axis=-1)

df = pd.DataFrame(f_out, columns=ai.names)

df['pos_x'] = grid.points[:,0]
df['pos_y'] = grid.points[:,1]

df = df.set_index(['pos_x', 'pos_y'])

df.to_csv('output/mdot0130_phi080_K100.csv')
# %%

df

# %%

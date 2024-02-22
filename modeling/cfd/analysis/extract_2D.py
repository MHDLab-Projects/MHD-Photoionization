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

sp_dir = gen_path('sharepoint')

results_dir = pjoin(sp_dir, 'Team Member Files', 'DaveH', 'Results', 'axiJP8200_17Jul23')
# results_dir = 'input'

fps = {
    '0.8_0.1': pjoin(results_dir, 'medium', 'mdot0130_phi080_K010', 'frontCyl_chem1.vtk'),
    '0.8_1.0': pjoin(results_dir, 'medium', 'mdot0130_phi080_K100', 'frontCyl_chem1.vtk'),
    '0.6_0.1': pjoin(results_dir, 'medium', 'mdot0130_phi060_K010', 'frontCyl_chem1.vtk'),
    '0.6_1.0': pjoin(results_dir, 'medium', 'mdot0130_phi060_K100', 'frontCyl_chem1.vtk'),
}



#%%

soi = ['K', 'Kp', 'em', 'OH', 'OHm', 'KOH', 'O2', 'H2O', 'N2', 'CO2']
soi_Yeq = ['Yeq_K', 'Yeq_K+', 'Yeq_e-', 'Yeq_OH', 'Yeq_OH-', 'Yeq_KOH', 'Yeq_K2CO3']
additional = ['T', 'p', 'rho']
all_fields = [*soi, *soi_Yeq, *additional]



# %%

from mhdpy.pyvista_utils import AxiMesh, pv_to_unstack_xr, downsel_arrays, pv_to_xr



#%%[markdown]

# Demo

#%%

# fp = pjoin(results_dir, 'medium', 'mdot0130_phi080_K100', 'frontCyl_chem1.vtk')
# fp = pjoin(results_dir, 'medium', 'mdot0130_phi080_K010', 'frontCyl_chem1.vtk')
fp = fps['0.8_1.0']

mesh = pv.read(fp)

# %%

mesh = downsel_arrays(mesh, [*soi_Yeq, *additional])

ds = pv_to_xr(mesh)['point']

df = ds.drop('pos_z').drop('dir').to_dataframe()

df

# pv_to_unstack_xr(mesh)

# %%



df.to_csv('output/mdot0130_phi080_K100.csv')
# %%

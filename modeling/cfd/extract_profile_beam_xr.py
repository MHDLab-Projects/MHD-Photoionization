#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from pvxarray.io import pyvista_to_xarray

from mhdpy.io import gen_path

from mhdpy.analysis.standard_import import *

pv.OFF_SCREEN = True
# %%

fname = 'output/torch_rot_interp.vtk'

m_torch = pv.read(fname)


m_torch.point_data.remove('U')
m_torch.point_data.remove('z')
m_torch.point_data.remove('y')

m_torch.slice(normal=(0,0,1)).plot()

# ds = pyvista_to_xarray(m_torch) # weird rastering problem...


#%%
from pv_utils import pv_to_xr

ds = pv_to_xr(m_torch)['point']
ds = ds.drop('dir')

ds2 = ds.set_index(i_point=('pos_x', 'pos_y','pos_z'), append=True)
ds2 = ds2.reset_index('i_point_level_0', drop=True).unstack('i_point')

ds2

#%%

ds2.sel(pos_z=0, method='nearest')['K'].plot()
# %%

ds2.sel(pos_x=.3, pos_y=0, method='nearest')['K'].plot()
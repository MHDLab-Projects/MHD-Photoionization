#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from mhdpy.io import gen_path

from mhdpy.analysis.standard_import import *

pv.OFF_SCREEN = True


r_beam = 0.005

# TODO: need to rotate frame or something in order to eventually convert to xarray. 
# zx_ratio = 1.875/4.25
# direc_beam = (-zx_ratio, 0, 1)

direc_beam = (0,0,1)

# %%

fname = 'output/torch_rot_interp.vtk'

m_torch = pv.read(fname)

m_torch.set_active_scalars('K')

m_torch

#%%

m_torch.slice_orthogonal().plot()

#%%

# https://docs.pyvista.org/version/stable/examples/01-filter/clipping.html

roi = pv.Cylinder(
    center=(0.25,0,0),
    direction=direc_beam,
    height=0.15,
    radius=r_beam
)

pl = pv.Plotter()
pl.add_mesh(roi, color='red')
pl.add_mesh(m_torch, opacity=0.5)

pl.show(cpos=('xy'))

#%%

# We only need to consider y direction within width of beam. 

m_torch_red = m_torch.clip_box([
   -np.inf, np.inf,
   -r_beam*1.1, r_beam*1.1,
   -np.inf, np.inf

], invert=False, crinkle=True)

#%%

pl = pv.Plotter()
pl.add_mesh(roi, color='red')
pl.add_mesh(m_torch_red, opacity=0.5)

pl.show(cpos=('xz'))

#%%

extracted = m_torch_red.clip_surface(roi, invert=True, crinkle=True)

#%%

# extracted.slice(
#     origin= (0, 0,0),
#     normal = (0,0,1)
# ).plot()

extracted.plot()


#%%

from mhdpy.pyvista import pv_to_unstack_xr

ds = pv_to_unstack_xr(extracted)



#%%

beam_xs = np.linspace(m_torch_red.bounds[0],m_torch_red.bounds[1], 20)


# TODO: having to interpolate on a common grid later due to rounding errors

ds_full = pv_to_unstack_xr(m_torch_red)
xs_full = ds_full.coords['pos_x'].values 
zero_ref_xs = xs_full - xs_full[0]

zero_ref_xs

#%%


dss = []
for x_beam in beam_xs:

    print("Beam Position: {}".format(x_beam))

    roi = pv.Cylinder(
        center=(x_beam,0,0),
        direction=direc_beam,
        height=0.15,
        radius=r_beam
    )

    extracted = m_torch_red.clip_surface(roi, invert=True, crinkle=True)

    ds = pv_to_unstack_xr(extracted)

    ds = ds.assign_coords(x_beam=x_beam)

    pos_xs = ds.coords['pos_x'].values 

    ds = ds.assign_coords(pos_x = pos_xs- pos_xs[0])


    ds = ds.interp(pos_x=zero_ref_xs, method='nearest')

    dss.append(ds)

#%%

ds = xr.concat(dss, 'x_beam')

ds = ds.dropna('pos_x','all')

# ds.mean('pos_z')['K'].plot()
# %%


g = ds['K'].mean('pos_z').plot(col='x_beam', col_wrap=4, sharey=False)




#%%

# dss[5]['K'].mean('pos_z')#.dropna('pos_x','all')

beam_conv = ds.mean('pos_z').mean('pos_x').mean('pos_y')

beam_conv.to_netcdf('output/beam_conv.cdf')


#%%

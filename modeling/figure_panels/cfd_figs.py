#%%

from mhdpy.analysis.standard_import import *

import pyvista as pv

pv.OFF_SCREEN = True


# %%

folder_cfd = pjoin(REPO_DIR, 'modeling','cfd')

fname = pjoin(folder_cfd, 'output','torch_rot_interp.vtk')

m_torch = pv.read(fname)

m_torch.set_active_scalars('K')

m_torch


# %%

ms = m_torch.slice(normal='y')

p = pv.Plotter()

p.camera_position = (0,1,0)

p.add_mesh(ms)

p.show()

#%%

from mhdpy.pyvista import pv_to_unstack_xr


ds = pv_to_unstack_xr(ms)
ds = ds.isel(pos_y=0)

# %%
ds['K'].plot()

#%%

dx = ds.coords['pos_x'].max() - ds.coords['pos_x'].min()
dx = dx.item()

dz = ds.coords['pos_z'].max() - ds.coords['pos_z'].min()
dz = dz.item()

aspect = dx/dz
aspect

#%%

output_dir = pjoin(DIR_FIG_OUT, 'cfd')
if not os.path.exists(output_dir): os.makedirs(output_dir)

for species in ['CO2','H2O','K','KOH', 'T']:
    # plt.figure(figsize= (2, 2*aspect))
    plt.figure()

    ds[species].plot()

    ax = plt.gca()
    ax.set_aspect(aspect)
    ax.set_title('')
    plt.tight_layout()
    plt.savefig(pjoin(output_dir,'{}.png'.format(species)))

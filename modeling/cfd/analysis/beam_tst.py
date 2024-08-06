#%%
from mhdpy.fileio import gen_path

from mhdpy.analysis.standard_import import *
import pint_pandas

import pyvista as pv
from mhdpy.pyvista_utils import AxiMesh


sys.path.append(pjoin(REPO_DIR, 'modeling', 'cfd'))
from pv_axi_utils import AxiInterpolator

sp_dir = gen_path('sharepoint')

from pi_paper_utils.fileio import cfd_all_fields, cfd_fp_dict


# %%

mesh = pv.read(cfd_fp_dict['0.8_1.0'])


#%%

import tri_mesh
import beam_utils

x_exit = 207.910 * 1e-3
x_min = 0.25
x_center = x_min
offset = np.array([-0.05,0.0,-0.1])

pos_center = np.array([x_center,0,0])
pos_source = pos_center + offset
pos_target = pos_center - offset

beam_axis = -offset
jet_axis = (1,0,0)
shat = beam_utils.coordinate_tensor(beam_axis, jet_axis)
L_beam = beam_utils.mag(offset)*2.0

x_beam = np.linspace(0,L_beam)
d_beam = 1e-2
# creates pos[:,:], mesh:pv.mesh, ds:xr.dataset, 

tri_out = tri_mesh.new_cylinder(x=x_beam,r=d_beam/2.0,n_theta=60, angle=00)

tri_out

#%%



#%%

ds = tri_out['ds']

def transform_beam(ds, shat):
    ds["pos_beam"] = ds["pos"]*1.0
    ds["pos"][:] = np.matmul(ds["pos_beam"].to_numpy(),shat) + pos_source

transform_beam(ds, shat)

#%%
# beam_utils.calc_intensity(ds)
mesh1 = beam_utils.extruded_mesh_from_dataset(ds, tri_out["tri"], add_bc=False)

bm = mesh1['beam']

#%%

p = pv.Plotter()

p.add_mesh(bm)
p.add_mesh(mesh)

p.show()

# %%
ai = AxiInterpolator(mesh, var_names = ['T','K'])

ai
#%%

f_out = ai(ds["pos"])
ds["T"] = ("s","ray"), f_out[:,:,0].T
ds["K"] = ("s","ray"), f_out[:,:,1].T

#%%

# beam_utils.calc_intensity(ds)
mesh1 = beam_utils.extruded_mesh_from_dataset(ds, tri_out["tri"], add_bc=False)

bm = mesh1['beam']

bm.set_active_scalars('T')

#%%

bm.plot()
# %%

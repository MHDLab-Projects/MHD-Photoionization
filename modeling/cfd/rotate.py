#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from mhdpy.io import gen_path

from mhdpy.analysis.standard_import import *

sp_dir = gen_path('sharepoint')

results_dir = pjoin(sp_dir, r"Team Member Files\DaveH\Results\axiJP8200_17Jul23")


#%%

case = "mdot0130_phi080_K010"
fname = os.path.join(results_dir, "medium",case,"frontCyl.vtk")

mesh = pv.read(fname)



for k in mesh.cell_data.keys():
    if k != 'em': mesh.cell_data.remove(k)


for k in mesh.point_data.keys():
    if k != 'em': mesh.point_data.remove(k)

mesh.set_active_scalars('em')

mesh


# %%


bounds = [0.21, 0.4,0,0.015, -0.1,0.1]
mesh_clip = mesh.clip_box(bounds, invert=False).extract_surface()

mesh_clip.plot(cpos='xy')


#%%

m3 = mesh_clip.extrude_rotate(
    angle=30,
    rotation_axis=(1,0,0)
    )


m3.plot(cpos='xy')


#%%

n_rots = 10

angles = np.linspace(0,360, n_rots)
angles = angles[:-1]

ms = []

for angle in angles:
    m = mesh_clip.rotate_x(angle)

    ms.append(m)



m_rot = pv.CompositeFilters.combine(ms)


m_rot.cell_data.remove('em')


m_rot = m_rot.cast_to_pointset()

m_rot.plot()

#%%


#%%


# x = np.linspace

bounds = m_rot.GetBounds()

spacing = 0.001

l_x = bounds[1] - bounds[0]
n_x = int(l_x/spacing)
l_y = bounds[3] - bounds[2]
n_y = int(l_y/spacing)
l_z = bounds[5] - bounds[4]
n_z = int(l_z/spacing)

# x = np.linspace(bounds[0],bounds[1], n_points)
# y = np.linspace(bounds[2],bounds[3], n_points)
# z = np.linspace(bounds[4],bounds[5], n_points)

# x, y, z = np.meshgrid(x,y,z)

# grid = pv.StructuredGrid(x,y,z)

# mi = grid.interpolate(m_rot)

grid = pv.ImageData()
grid.origin = (bounds[0], bounds[2], bounds[4])
grid.spacing = (spacing, spacing, spacing)
grid.dimensions = (n_x, n_y, n_z)

# grid = grid.cast_to_rectilinear_grid()

#%%

grid.slice_orthogonal().plot(show_vertices=True)

#%%

grid.plot(show_vertices=True, opacity=0.1)

# %%

p = pv.Plotter()

p.add_mesh(grid, show_vertices=True, opacity=0.1)
p.add_mesh(m_rot)

p.show()

#%%

m_i = grid.interpolate(m_rot, radius=spacing)


# m_i.plot()

#%%
p = pv.Plotter()

p.add_mesh(m_i, show_vertices=True, opacity=0.1)
# p.add_mesh(m_rot, show_vertices=True, opacity=0.1)

p.show()



#%%

m_i.slice_orthogonal().plot()
# %%

# m_i.save('output/em_mesh.vtk')

# m_out = m_i.cast_to_poly_points()

# m_out.save('output/em_mesh.ply', binary=False)

# m_i.save('output/test.csv', binary=False)

m_i.save('output/test.vtk')

# pv.save_meshio(m_out, 'output/em_mesh.stl')

#%%

import pandas as pd

data = dict(zip(['x','y','z'], m_i.points.T))

df = pd.DataFrame(data)

df['em'] = m_i['em']

df['em'] = df['em']/df['em'].max() # Temporary

df.to_csv('output/em_mesh.csv', index=False, header=False)



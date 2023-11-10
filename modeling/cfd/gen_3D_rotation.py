#%%
"""
Script to create a 3D ImageData vtk file from a 2D axi-symmetric simulation. The
2D plane is rotated a number of times, then interpolated to a 3D grid. The final
spacing of the grid is used to determine the needed angle of rotation and
interpolation radius. The data of a specified array is also output as a csv file
suitable for loading into comsol
"""

#TODO: the algorithm does not seem to be working well at long lengths from the torch. Could be related to a lower mesh size there. 
#TODO: The interpolation at the edges of the square, outside the rotated circle, is not defined. Perhaps reduce dimensions of imagedata mesh to be inside the rotated cylinder

#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os
import pandas as pd

from mhdpy.fileio import gen_path

from mhdpy.analysis.standard_import import *

pv.OFF_SCREEN = True

#%%

sp_dir = gen_path('sharepoint')
results_dir = pjoin(sp_dir, r"Team Member Files\DaveH\Results\axiJP8200_17Jul23")
case = "mdot0130_phi080_K010"
fname = os.path.join(results_dir, "medium",case,"frontCyl.vtk")

mesh = pv.read(fname)



mesh.set_active_scalars('CO2')

mesh


# %%

# Clip the dataset to this box. 
bounds = [0.21, 0.5,0,0.05, -0.1,0.1]

l_x = bounds[1] - bounds[0]
l_y = bounds[3] - bounds[2]
l_z = bounds[5] - bounds[4]

mesh_clip = mesh.clip_box(bounds, invert=False).extract_surface()

mesh_clip.plot(cpos='xy')

#%%

# Ultimate spacing of the 3D grid. This will also be used for the interpolation radius and rotation angle determineation
spacing = 0.001


# We need the rotational angle to fine enough that the spacing at the end of the arc is equal to this
# Tan(theta) = spacing/L ~ theta
# n_rots = 2pi/theta = 2*pi*L/spacing

L = l_y

n_rots = int(np.ceil(2*np.pi*L/spacing))

print("Generating {} rotations".format(n_rots))

angles = np.linspace(0,360, n_rots)
angles = angles[:-1]

ms = []

for angle in angles:
    m = mesh_clip.rotate_x(angle)

    ms.append(m)

m_rot = pv.CompositeFilters.combine(ms)

# m_rot.cell_data.remove('em')

m_rot = m_rot.cast_to_pointset()

#%%
m_rot.plot()

#%%

bounds = m_rot.GetBounds()

l_x = bounds[1] - bounds[0]
l_y = bounds[3] - bounds[2]
l_z = bounds[5] - bounds[4]

n_x = int(l_x/spacing)
n_y = int(l_y/spacing)
n_z = int(l_z/spacing)

grid = pv.ImageData()
grid.origin = (bounds[0], bounds[2], bounds[4])
grid.spacing = (spacing, spacing, spacing)
grid.dimensions = (n_x, n_y, n_z)

# grid = grid.cast_to_rectilinear_grid()

# %%

p = pv.Plotter()

p.add_mesh(grid, show_vertices=True, opacity=0.1)
p.add_mesh(m_rot)

p.show()

#%%

m_i = grid.interpolate(m_rot, radius=spacing)

m_i.plot()


#%%

# m_i.slice_orthogonal().plot()
# %%

m_i.save('output/torch_rot_interp.vtk')

#%%


data = dict(zip(['x','y','z'], m_i.points.T))

df = pd.DataFrame(data)

df['em'] = m_i['em']

df['em'] = df['em']/df['em'].max() # Temporary

df.to_csv('output/em_mesh.csv', index=False, header=False)



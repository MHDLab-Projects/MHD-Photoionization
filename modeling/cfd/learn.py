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

mesh


#%%

points = np.random.rand(100,3)

mesh = pv.PolyData(points)

mesh.plot()

#%%



mesh = pv.examples.load_hexbeam()

mesh.plot(show_edges=True)

#%%
mesh.active_scalars_name = 'sample_point_scalars'

mesh.plot(show_edges=True)

#%%

np_points = np.array([[0, 0, 0],
                      [1, 0, 0],
                      [0.5, 0.667, 0]])

cells = [3, 0, 1, 2]

mesh = pv.PolyData(np_points, cells)

mesh.plot(cpos='xy', point_size=20)

#%%

grid = pv.ImageData(dimensions=(3,3,1))

ugrid = grid.cast_to_unstructured_grid()

ugrid.plot(show_edges=True)

#%%

simple_range = range(ugrid.n_cells)

ugrid.cell_data['cdata'] = simple_range

ugrid.plot()

#%%

ugrid.point_data['pdata'] = range(ugrid.n_points)

ugrid.point_data
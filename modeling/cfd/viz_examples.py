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
# print("array names")
# print(mesh.array_names)


#%%
p = pv.Plotter()
# mesh['K']
p.add_mesh(mesh, scalars='K', clim=[0,0.0001])

p.camera_position = 'xy'
p.show()

#%%
cell = mesh.cell_centers()
pos= cell.points
x = cell["CO"]
y = (x - x.min())/(x.max() - x.min())
plt.pcolor(pos[:,0], pos[:,1], x)

#%%
unit_cm = 1e-2
x_exit = 20.8*unit_cm

var_names = "OH T K KOH Kp em".split()
locs = [-1, 0, 1, 5, 10, 20, 40]
n_locs = len(locs)
n_vars = len(var_names)
fig, ax = plt.subplots(n_vars, n_locs, sharey=True)
for i_loc, x_pos_cm in enumerate(locs):
    a = [x_exit + x_pos_cm*unit_cm, 0,    0]
    b = [x_exit + x_pos_cm*unit_cm, 2*unit_cm, 0]
    line1 = mesh.sample_over_line(a, b)    
    for i_var, var_name in enumerate(var_names):    
        ax[i_var, i_loc].plot(line1[var_name][:], line1.points[:,1]/unit_cm)
        
        ax[0, i_loc].set_title("x - x_exit = {} [cm]".format(x_pos_cm))
        ax[i_var, i_loc].set_xlabel("{}".format(var_name))
        
fig.set_size_inches(12,8)
fig.tight_layout()

#%%
locs = [1e-3, 5e-3]
n_vars = len(var_names)
n_locs = len(locs)
fig, ax = plt.subplots(n_vars, n_locs, sharex=True)
for i_loc, y_pos_cm in enumerate(locs):
    a = [x_exit - 1e-2 , y_pos_cm*unit_cm , 0]
    b = [x_exit + 40e-2, y_pos_cm*unit_cm, 0]
    line1 = mesh.sample_over_line(a, b)    
    for i_var, var_name in enumerate(var_names):    
        ax[i_var, i_loc].plot( (line1.points[:,0] - x_exit)/unit_cm, line1[var_name][:])
        
        ax[0, i_loc].set_title("y = {} [cm]".format(y_pos_cm))
        
        ax[i_var, 0].set_ylabel(var_name)
        ax[-1, i_loc].set_xlabel("{}".format("x - x_exit [cm]"))

fig.set_size_inches(12,8)
fig.tight_layout()

#%%

mesh['K']

#%%

# from pyvista import examples

# mesh = examples.load_airplane()
# # mesh.plot()

# mesh
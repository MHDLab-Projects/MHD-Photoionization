#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from mhdpy.io import gen_path

from mhdpy.analysis.standard_import import *

sp_dir = gen_path('sharepoint')

results_dir = pjoin(sp_dir, r"Team Member Files\DaveH\Results\axiJP8200_17Jul23")
# %%
case = "mdot0130_phi080_K010"
fname = os.path.join(results_dir, "medium",case,"frontCyl.vtk")

mesh = pv.read(fname)

mesh
# %%
unit_cm = 1e-2
x_exit = 20.8*unit_cm

y_pos_cm = 0

a = [x_exit - 1e-2 , y_pos_cm*unit_cm , 0]
b = [x_exit + 40e-2, y_pos_cm*unit_cm, 0]
line1 = mesh.sample_over_line(a, b)    

#%%

soi = ['K', 'Kp', 'em', 'OH', 'OHm', 'KOH']

for species in soi:

    plt.plot(line1.points[:,0], line1[species], label=species)

plt.yscale('log')
plt.ylabel("Species concentration")
plt.xlabel("X Position")
plt.legend()

#%%

data = np.array([
    line1[species] for species in soi
])

data = data.T

df_out = pd.DataFrame(data, index = line1.points[:,0], columns=soi)
df_out.index.name = 'x'

df_out.to_csv('output/line_profies.csv')
# %%

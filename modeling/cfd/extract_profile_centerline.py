#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from mhdpy.fileio import gen_path

from mhdpy.analysis.standard_import import *

sp_dir = gen_path('sharepoint')

results_dir = pjoin(sp_dir, 'Team Member Files', 'DaveH', 'Results', 'axiJP8200_17Jul23')

# %%
def extract_line(mesh, fields):

    unit_cm = 1e-2
    x_exit = 20.8*unit_cm

    y_pos_cm = 0

    a = [x_exit - 1e-2 , y_pos_cm*unit_cm , 0]
    b = [x_exit + 40e-2, y_pos_cm*unit_cm, 0]
    line1 = mesh.sample_over_line(a, b)    

    data = np.array([
        line1[species] for species in fields
    ])

    data = data.T

    df_out = pd.DataFrame(data, index = line1.points[:,0], columns=all_fields)
    df_out.index.name = 'x'

    return df_out

soi = ['K', 'Kp', 'em', 'OH', 'OHm', 'KOH']
additional = ['T', 'p']
all_fields = [*soi, *additional]

#%%


cases = [
    "mdot0130_phi080_K000",
    "mdot0130_phi080_K010",
    "mdot0130_phi080_K100"
]

kwts = [0,0.1,1] #TODO: regex


dfs = []
for i, case in enumerate(cases):
    fname = os.path.join(results_dir, "medium",case,"frontCyl.vtk")

    mesh = pv.read(fname)

    df_out = extract_line(mesh, all_fields)

    df_out['kwt'] = kwts[i]

    dfs.append(df_out)


df_lines = pd.concat(dfs)
df_lines = df_lines.set_index('kwt', append=True)
#%%

ds = xr.Dataset.from_dataframe(df_lines)

ds

#%%

ds[soi].to_array('species').plot(row='kwt', hue='species')

plt.yscale('log')
plt.ylabel("Species concentration")
plt.xlabel("X Position")


#%%



ds.to_netcdf('output/line_profiles.cdf')
# %%

#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from mhdpy.fileio import gen_path

from mhdpy.analysis.standard_import import *

sp_dir = gen_path('sharepoint')

# results_dir = pjoin(sp_dir, 'Team Member Files', 'DaveH', 'Results', 'axiJP8200_17Jul23')
results_dir = 'input'

# %%
def extract_line(mesh, fields, a, b):

    line1 = mesh.sample_over_line(a, b)    

    data = np.array([
        line1[species] for species in fields
    ])

    data = data.T

    df_out = pd.DataFrame(data, columns=all_fields)
    df_out['x'] = line1.points[:,0]
    df_out['y'] = line1.points[:,1]
    df_out['z'] = line1.points[:,2]

    df_out = df_out.set_index(['x', 'y', 'z'])
    # df_out.index.name = 'x'

    return df_out

soi = ['K', 'Kp', 'em', 'OH', 'OHm', 'KOH']
additional = ['T', 'p']
all_fields = [*soi, *additional]

#%%

unit_cm = 1e-2
x_exit = 20.8*unit_cm

y_pos_cm = 0

a = [x_exit - 1e-2 , y_pos_cm*unit_cm , 0]
b = [x_exit + 40e-2, y_pos_cm*unit_cm, 0]

cases = [
    "mdot0130_phi080_K000",
    "mdot0130_phi080_K010",
    "mdot0130_phi080_K100"
]

kwts = [0,0.1,1] #TODO: regex


dfs = []
for i, case in enumerate(cases):
    fname = os.path.join(results_dir, "medium",case,"frontCyl.vtk")
    # fname = os.path.join(results_dir, case,"frontCyl.vtk")

    mesh = pv.read(fname)



    df_out = extract_line(mesh, all_fields, a, b)
    df_out = df_out.reset_index('y').reset_index('z').drop('y', axis=1).drop('z', axis=1)   

    df_out['kwt'] = kwts[i]

    dfs.append(df_out)


df_lines = pd.concat(dfs)
df_lines = df_lines.set_index('kwt', append=True)

#%%

df_lines
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
# %%[markdown]

# # Extracting profiles along beam path


def gen_beam_line(position):
    xy_ratio = 1.875/4.25

    unit_cm = 1e-2
    x_exit = 20.8*unit_cm

    #position of interseciton of line with z=0
    x_beam = x_exit +  position*unit_cm

    # make two points that intersect this point and are in the direction of the beam
    y_distance = 5*unit_cm

    a = [x_beam , 0 , 0]
    b = [x_beam + y_distance*xy_ratio, y_distance, 0]

    return a, b

a,b = gen_beam_line(30)

line = pv.Line(a, b, resolution=100)

# pv.plot([line, mesh])

#%%

tc = '536_pos'
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

beam_positions = ds_absem.coords['motor'].pint.quantify('mm')

beam_positions = beam_positions.pint.to('cm').pint.magnitude

beam_positions

#%%




dfs = []

for pos in beam_positions:

    a, b = gen_beam_line(pos)

    df = extract_line(mesh, all_fields, a, b)
    df = df.reset_index('z').drop('z', axis=1)


    df = df.reset_index('x').reset_index('y')
    df['x'] = df['x'] - df['x'].min()
    df['y'] = df['y'] - df['y'].min()

    df['dist'] = np.sqrt(df['x']**2 + df['y']**2)

    df['pos'] = pos

    df = df.set_index(['pos', 'dist'])

    dfs.append(df)
    # df = df.set_index('dist')

df_beam = pd.concat(dfs)
#%%

ds = xr.Dataset.from_dataframe(df_beam)

ds

#%%

from mhdpy.plot import dropna
g = ds['K'].plot(hue='pos')


plt.yscale('log')

plt.ylim(1e-8, 1e-2)

dropna(g)

# df['K'].plot()

#%%

ds['K'].sel(pos=25, method='nearest').dropna('dist').plot()



#%%

ds['K'].sel(dist=0).plot()
plt.yscale('log')

#%%

ds_out = ds[['K', 'T', 'p']]



ds_out.to_netcdf('output/nK_beam_profiles.cdf')





# %%

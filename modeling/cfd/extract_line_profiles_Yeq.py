"""
Script to extract line and beam profiles from CFD simulation data with new Yeq fields
TODO: replace old extract_line_profiles.py with this script. Waiting on K=1% data. 
"""


#%%
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import os

from mhdpy.fileio import gen_path

from mhdpy.analysis.standard_import import *

sp_dir = gen_path('sharepoint')

results_dir = pjoin(sp_dir, 'Team Member Files', 'DaveH', 'Results', 'axiJP8200_17Jul23')
# results_dir = 'input'

fps = {
    '0.1': pjoin(results_dir, 'medium', 'mdot0130_phi080_K010', 'frontCyl_chem1.vtk'),
    '1.0': pjoin(results_dir, 'medium', 'mdot0130_phi080_K100', 'frontCyl_chem1.vtk')
}




#%%

soi = ['K', 'Kp', 'em', 'OH', 'OHm', 'KOH', 'O2', 'H2O', 'N2', 'CO2']
soi_Yeq = ['Yeq_K', 'Yeq_K+', 'Yeq_e-', 'Yeq_OH', 'Yeq_OH-', 'Yeq_KOH', 'Yeq_K2CO3']
additional = ['T', 'p']
all_fields = [*soi, *soi_Yeq, *additional]



# %%

from pv_axi_utils import AxiInterpolator,AxiMesh


def gen_beam_line(position, xy_ratio=1.875/4.25):
    # xy_ratio=1

    unit_cm = 1e-2
    x_exit = 20.8*unit_cm

    #position of interseciton of line with z=0
    x_beam = x_exit +  position*unit_cm

    # make two points that intersect this point and are in the direction of the beam
    y_distance = 5*unit_cm

    a = [x_beam - y_distance*xy_ratio, y_distance , 0]
    b = [x_beam + y_distance*xy_ratio, -y_distance, 0]

    return a, b


def extract_line_axi(mesh, a, b):

    am = AxiMesh(mesh)

    line1 = am.sample_over_line(a, b, resolution=100)

    line_cart  = pv.Line(a, b, 100)

    line1.points = line_cart.points

    # calculate the distance and values along the line

    return line1

def convert_line_df(line, fields):
    dist = np.zeros(line.n_points)

    for i in range(1, line.n_points):
        dist[i] = dist[i-1] + np.linalg.norm(line.points[i] - line.points[i-1])

    lines_out = []
    for field in fields:
        vals = line.point_data[field]

        line_out = pd.Series(vals, index=dist)

        lines_out.append(line_out)

    df = pd.concat(lines_out, axis=1)
    df.columns = fields
    return df

def interp_df_to_new_index(df_out, new_index):
    # Initialize a new DataFrame with the new index
    df_interpolated = pd.DataFrame(index=new_index)

    # Interpolate each column separately
    for column in df_out.columns:
        df_interpolated[column] = np.interp(new_index, df_out.index, df_out[column])

    return df_interpolated

#%%[markdown]

# Demo

#%%

# fp = pjoin(results_dir, 'medium', 'mdot0130_phi080_K100', 'frontCyl_chem1.vtk')
# fp = pjoin(results_dir, 'medium', 'mdot0130_phi080_K010', 'frontCyl_chem1.vtk')
fp = fps['1.0']

mesh = pv.read(fp)

a, b = gen_beam_line(3, xy_ratio=1.875/4.25)

line1 = extract_line_axi(mesh, a, b)

p = pv.Plotter()

p.add_mesh(mesh, scalars='T')
p.add_mesh(line1, color='red')

p.camera_position = [(0, 0, 1), (0.1, 0, 0), (0, 0, 0)]

p.show()

#%%

df_lines = convert_line_df(line1, ['Yeq_K'])

df_lines.plot()

# plt.yscale('log')

midpoint = max(df_lines.index)/2
plt.axvline(midpoint, color='gray', linestyle='--')

#%%[markdown]

# # Extract lines along beam axis


#%%
tc = '536_pos'
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

beam_positions = ds_absem.coords['motor'].pint.quantify('mm')

beam_positions = beam_positions.pint.to('cm').pint.magnitude

beam_positions


# beam_positions = np.arange(1, 30, 5)
dist_grid = np.arange(0, 0.1, 0.001)

dss = []


for kwt, fp in fps.items():

    mesh = pv.read(fp)
    for position in beam_positions:

        a, b = gen_beam_line(position)

        line_out = extract_line_axi(mesh, a, b)

        df_lines = convert_line_df(line_out, ['T','p','Yeq_K'])

        df_int = interp_df_to_new_index(df_lines, dist_grid)

        ds_out = xr.Dataset(df_int)
        ds_out = ds_out.rename({'dim_0':'dist'})

        ds_out = ds_out.assign_coords(pos=position)
        ds_out = ds_out.assign_coords(kwt=float(kwt))

        ds_out.expand_dims('pos').expand_dims('kwt').stack(temp=('pos', 'kwt'))

        dss.append(ds_out)

    # break

ds_lines = xr.concat(dss, dim='temp')

ds_lines = ds_lines.set_index(temp=['pos', 'kwt']).unstack('temp')

#%%

ds_lines
#%%

ds_lines['T'].plot(col='kwt', hue='pos')

#%%

ds_lines['Yeq_K'].plot(col='kwt', hue='pos')

plt.yscale('log')

plt.ylim(1e-14, 1e-2)

#%%

ds_lines.to_netcdf(pjoin('output', 'line_profiles_beam_Yeq.cdf'))


#%%[markdown]

# # Extract lines along torch axis

#%%


unit_cm = 1e-2
x_exit = 20.8*unit_cm

a = [x_exit - 1e-2 , 0 , 0]
b = [x_exit + 40e-2, 0, 0]

line1 = extract_line_axi(mesh, a, b)

p = pv.Plotter()

p.add_mesh(mesh, scalars='K')
p.add_mesh(line1, color='red', line_width=5)


p.camera_position = [(0, 0, 1), (0.1, 0, 0), (0, 0, 0)]


p.show()

#%%

#TODO: add multi file for beam positions above. 


dss = []

for kwt, fp in fps.items():

    mesh = pv.read(fp)

    dist_grid = np.arange(0, 0.4, 0.001)

    line_out = extract_line_axi(mesh, a, b)

    df_lines = convert_line_df(line_out, all_fields)

    df_int = interp_df_to_new_index(df_lines, dist_grid)

    ds_out = xr.Dataset(df_int)
    ds_out = ds_out.rename({'dim_0':'x'})

    ds_out
    ds_out = ds_out.assign_coords(kwt=float(kwt))

    dss.append(ds_out)

    # ds_out = ds_out.drop_vars(['K', 'Kp'])
    # ds_out = ds_out.rename({'Yeq_K':'K', 'Yeq_K+':'Kp'})

    # # 

ds_out = xr.concat(dss, 'kwt')

ds_out.to_netcdf(pjoin('output', 'line_profiles_torchaxis_Yeq.cdf'))


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
import pint_pandas

sp_dir = gen_path('sharepoint')

results_dir = pjoin(sp_dir, 'Team Member Files', 'DaveH', 'Results', 'axiJP8200_17Jul23')
# results_dir = 'input'

fps = {
    '0.8_0.1': pjoin(results_dir, 'medium', 'mdot0130_phi080_K010', 'frontCyl_chem1.vtk'),
    '0.8_1.0': pjoin(results_dir, 'medium', 'mdot0130_phi080_K100', 'frontCyl_chem1.vtk'),
    '0.6_0.1': pjoin(results_dir, 'medium', 'mdot0130_phi060_K010', 'frontCyl_chem1.vtk'),
    '0.6_1.0': pjoin(results_dir, 'medium', 'mdot0130_phi060_K100', 'frontCyl_chem1.vtk'),
}



#%%
tc = '536_pos'
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

beam_positions = ds_absem.coords['motor'].pint.quantify('mm')

beam_positions = beam_positions.pint.to('cm').pint.magnitude
beam_positions = [Quantity(pos, 'cm') for pos in beam_positions]

beam_positions

#%%

#offsets in cm (TODO: pint)
beam_offsets = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
beam_offsets = [Quantity(offset, 'cm') for offset in beam_offsets]



#%%

soi = ['K', 'Kp', 'em', 'OH', 'OHm', 'KOH', 'O2', 'H2O', 'N2', 'CO2']
soi_Yeq = ['Yeq_K', 'Yeq_K+', 'Yeq_e-', 'Yeq_OH', 'Yeq_OH-', 'Yeq_KOH', 'Yeq_K2CO3']
additional = ['T', 'p', 'rho']
all_fields = [*soi, *soi_Yeq, *additional]



# %%

from mhdpy.pyvista_utils import AxiMesh


def gen_beam_line(
        ta_position: Quantity, 
        ta_offset: Quantity, 
        xy_ratio=Quantity(1.875/4.25, 'cm/cm'),
        ):
    """
    Generate a line with respect to torch axis
    ta_position: position down the torch axis
    ta_offset: offset from the torch axis (z) 
    """

    x_exit = Quantity(20.8, 'cm')

    #position of interseciton of line with z=0
    x_beam = x_exit +  ta_position

    # make two points that intersect this point and are in the direction of the beam
    y_distance = Quantity(5, 'cm')

    a = [x_beam - y_distance*xy_ratio, y_distance , ta_offset]
    b = [x_beam + y_distance*xy_ratio, -y_distance, ta_offset]

    return a, b


def extract_line_axi(mesh, a: Quantity, b:Quantity):

    am = AxiMesh(mesh)


    a = [e.to('m').magnitude for e in a]
    b = [e.to('m').magnitude for e in b]

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
fp = fps['0.8_1.0']

mesh = pv.read(fp)

a, b = gen_beam_line(Quantity(3, 'cm'), Quantity(0,'cm'), xy_ratio=Quantity(1.875/4.25, 'cm/cm'))

line1 = extract_line_axi(mesh, a, b)

p = pv.Plotter()

p.add_mesh(mesh, scalars='T')
p.add_mesh(line1, color='red')

p.camera_position = [(0, 0, 1), (0.1, 0, 0), (0, 0, 0)]

# p.show()

#%%

df_lines = convert_line_df(line1, ['Yeq_K'])

df_lines.plot()

# plt.yscale('log')

midpoint = max(df_lines.index)/2
plt.axvline(midpoint, color='gray', linestyle='--')

#%%[markdown]

# # Extract lines along beam axis

# beam_positions = np.arange(1, 30, 5)
dist_grid = np.arange(0, 0.1, 0.001)
dist_grid = Quantity(dist_grid, 'm')

dss = []


for key, fp in fps.items():

    phi, kwt = key.split('_')

    mesh = pv.read(fp)
    for position in beam_positions:
        for offset in beam_offsets:

            a, b = gen_beam_line(position, offset)

            line_out = extract_line_axi(mesh, a, b)

            #TODO: extract all soi
            df_lines = convert_line_df(line_out, ['T','p','rho', 'Yeq_K'])

            df_int = interp_df_to_new_index(df_lines, dist_grid.to('m').magnitude)

            ds_out = xr.Dataset(df_int)
            ds_out = ds_out.rename({'dim_0':'dist'})

            ds_out = ds_out.assign_coords(pos=position.to('m').magnitude)
            ds_out = ds_out.assign_coords(kwt=float(kwt))
            ds_out = ds_out.assign_coords(phi=float(phi))
            ds_out = ds_out.assign_coords(offset=offset.to('m').magnitude)

            ds_out.expand_dims('pos').expand_dims('kwt').expand_dims('phi').expand_dims('offset').stack(temp=('pos', 'kwt', 'phi', 'offset'))

            dss.append(ds_out)

    # break

ds_lines = xr.concat(dss, dim='temp')

ds_lines = ds_lines.set_index(temp=['pos', 'kwt', 'phi', 'offset']).unstack('temp')

#%%

ds_lines['pos'] = ds_lines['pos'].pint.quantify('m').pint.to('cm')
ds_lines['offset'] = ds_lines['offset'].pint.quantify('m').pint.to('cm')
ds_lines['dist'] = ds_lines['dist'].pint.quantify('m').pint.to('m')

#%%

ds_lines

#%%

ds_lines['T'].sel(offset=0).plot(col='kwt', hue='pos', row='phi')

#%%


ds_lines['T'].sel(kwt=1).plot(col='phi', hue='pos', row='offset')


#%%

ds_lines['Yeq_K'].sel(offset=0).plot(col='kwt', hue='pos', row='phi')

plt.yscale('log')

plt.ylim(1e-14, 1e-2)

#%%

ds_lines['Yeq_K'].sel(kwt=1).plot(col='phi', hue='pos', row='offset')
plt.yscale('log')


#%%

ds_lines.pint.dequantify().to_netcdf(pjoin('output', 'line_profiles_beam_Yeq.cdf'))


#%%[markdown]

# # Extract lines along torch axis

#%%


x_exit = Quantity(20.8, 'cm')

a = [x_exit - Quantity(1, 'cm') , Quantity(0, 'cm') , Quantity(0, 'cm')]
b = [x_exit + Quantity(40, 'cm'), Quantity(0, 'cm'), Quantity(0, 'cm')]

line1 = extract_line_axi(mesh, a, b)

p = pv.Plotter()

p.add_mesh(mesh, scalars='K')
p.add_mesh(line1, color='red', line_width=5)


p.camera_position = [(0, 0, 1), (0.1, 0, 0), (0, 0, 0)]


# p.show()

#%%

#TODO: add multi file for beam positions above. 


dss = []

for key, fp in fps.items():
    for offset in beam_offsets:

        #TODO: pint
        a = [x_exit - Quantity(1, 'cm') , Quantity(0, 'cm'), offset]
        b = [x_exit + Quantity(40, 'cm'), Quantity(0, 'cm'), offset]

        phi, kwt = key.split('_')

        mesh = pv.read(fp)

        dist_grid = np.arange(0, 0.4, 0.001)
        dist_grid = Quantity(dist_grid, 'm')

        line_out = extract_line_axi(mesh, a, b)

        df_lines = convert_line_df(line_out, all_fields)

        df_int = interp_df_to_new_index(df_lines, dist_grid.to('m').magnitude)

        ds_out = xr.Dataset(df_int)
        ds_out = ds_out.rename({'dim_0':'x'})

        ds_out
        ds_out = ds_out.assign_coords(kwt=float(kwt))
        ds_out = ds_out.assign_coords(phi=float(phi))
        ds_out = ds_out.assign_coords(offset=offset.to('m').magnitude)


        ds_out.expand_dims('kwt').expand_dims('phi').expand_dims('offset').stack(temp=('kwt', 'phi', 'offset'))

        dss.append(ds_out)

        # ds_out = ds_out.drop_vars(['K', 'Kp'])
        # ds_out = ds_out.rename({'Yeq_K':'K', 'Yeq_K+':'Kp'})

        # # 

ds_out = xr.concat(dss, 'temp')

ds_out = ds_out.set_index(temp=['kwt', 'phi', 'offset']).unstack('temp')

ds_out['x'] = ds_out['x'].pint.quantify('m').pint.to('m')
ds_out['offset'] = ds_out['offset'].pint.quantify('m').pint.to('cm')

ds_out.pint.dequantify().to_netcdf(pjoin('output', 'line_profiles_torchaxis_Yeq.cdf'))


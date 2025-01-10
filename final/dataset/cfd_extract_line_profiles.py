"""
Script to extract line and beam profiles from CFD simulation data with new Yeq fields
TODO: replace old extract_line_profiles.py with this script. Waiting on K=1% data. 
"""

#%%
from mhdpy.fileio import gen_path
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pint_pandas
import pyvista as pv
from mhdpy.pyvista_utils import AxiMesh
from pi_paper_utils.fileio import cfd_fp_dict, cfd_all_fields
import pi_paper_utils as ppu

from pi_paper_utils.constants import CFD_EXIT_OFFSET

BEAM_Y_DISTANCE = Quantity(5, 'cm')
beam_path_dist_grid = np.arange(0, 0.1, 0.0005)
beam_path_dist_grid = Quantity(beam_path_dist_grid, 'm')

#%%
tc = '536_pos'
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data')

ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem', '{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

beam_positions = ds_absem.coords['motor'].pint.quantify('mm')
beam_positions = beam_positions.pint.to('cm').pint.magnitude
beam_positions = [Quantity(pos, 'cm') for pos in beam_positions]

#%%
beam_offsets = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
beam_offsets = [Quantity(offset, 'cm') for offset in beam_offsets]

def gen_beam_line(ta_position: Quantity, ta_offset: Quantity, xy_ratio=Quantity(1.875/4.25, 'in/in')):
    """
    Generate a line with respect to torch axis
    ta_position: position down the torch axis
    ta_offset: offset from the torch axis (z)
    """
    x_beam = CFD_EXIT_OFFSET + ta_position
    a = [x_beam - BEAM_Y_DISTANCE * xy_ratio, BEAM_Y_DISTANCE, ta_offset]
    b = [x_beam + BEAM_Y_DISTANCE * xy_ratio, -BEAM_Y_DISTANCE, ta_offset]
    return a, b

def extract_line_axi(mesh, a: Quantity, b: Quantity, resolution=1000):
    """
    Create a Line object from a mesh object along the line defined by a and b
    """
    am = AxiMesh(mesh)
    a = [e.to('m').magnitude for e in a]
    b = [e.to('m').magnitude for e in b]
    line1 = am.sample_over_line(a, b, resolution=resolution)
    line_cart = pv.Line(a, b, resolution=resolution)
    line1.points = line_cart.points
    return line1

def convert_line_df(line, fields):
    """
    Converts a pyvista Line object into a dataframe with index distance along the line
    """
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
    df_interpolated = pd.DataFrame(index=new_index)
    for column in df_out.columns:
        df_interpolated[column] = np.interp(new_index, df_out.index, df_out[column])
    return df_interpolated

#%%[markdown]
# Demo

#%%
fp = cfd_fp_dict['0.8_1.0']
mesh = pv.read(fp)
a, b = gen_beam_line(Quantity(3, 'cm'), Quantity(0, 'cm'), xy_ratio=Quantity(1.875/4.25, 'cm/cm'))
line1 = extract_line_axi(mesh, a, b)

p = pv.Plotter()
p.add_mesh(mesh, scalars='T')
p.add_mesh(line1, color='red')
p.camera_position = [(0, 0, 1), (0.1, 0, 0), (0, 0, 0)]
# p.show()

#%%
df_lines = convert_line_df(line1, [ppu.CFD_K_SPECIES_NAME])
df_lines.plot()
midpoint = max(df_lines.index) / 2
plt.axvline(midpoint, color='gray', linestyle='--')

#%%[markdown]
# Extract lines along beam axis for different beam positions and offsets


dss = []

for key, fp in cfd_fp_dict.items():
    phi, kwt = key.split('_')
    mesh = pv.read(fp)
    for position in beam_positions:
        for offset in beam_offsets:
            a, b = gen_beam_line(position, offset)
            line_out = extract_line_axi(mesh, a, b)
            df_lines = convert_line_df(line_out, cfd_all_fields)
            df_int = interp_df_to_new_index(df_lines, beam_path_dist_grid.to('m').magnitude)
            ds_out = xr.Dataset(df_int)
            ds_out = ds_out.rename({'dim_0': 'dist'})
            ds_out = ds_out.assign_coords(motor=position.to('m').magnitude)
            ds_out = ds_out.assign_coords(kwt=float(kwt))
            ds_out = ds_out.assign_coords(phi=float(phi))
            ds_out = ds_out.assign_coords(offset=offset.to('m').magnitude)
            ds_out.expand_dims('motor').expand_dims('kwt').expand_dims('phi').expand_dims('offset').stack(temp=('motor', 'kwt', 'phi', 'offset'))
            dss.append(ds_out)

ds_lines = xr.concat(dss, dim='temp')
ds_lines = ds_lines.set_index(temp=['motor', 'kwt', 'phi', 'offset']).unstack('temp')

#%%
ds_lines['motor'] = ds_lines['motor'].pint.quantify('m').pint.to('mm')
ds_lines['offset'] = ds_lines['offset'].pint.quantify('m').pint.to('mm')
ds_lines['dist'] = ds_lines['dist'].pint.quantify('m').pint.to('mm')

#%%
ds_lines[ppu.CFD_K_SPECIES_NAME].sel(phi=0.8, kwt=1).sel(motor=100, method='nearest').sel(offset=0).dropna('dist')

#%%
ds_lines.pint.dequantify().to_netcdf(pjoin('output', 'cfd_profiles_beam_mobile.cdf'))

#%%[markdown]
# Extract barrel exit profile

#%%
exit_offset = Quantity(5, 'mm')
a = [CFD_EXIT_OFFSET, Quantity(-5, 'cm'), Quantity(0, 'cm')]
b = [CFD_EXIT_OFFSET, Quantity(5, 'cm'), Quantity(0, 'cm')]
line1 = extract_line_axi(mesh, a, b)

p = pv.Plotter()
p.add_mesh(mesh, scalars='K')
p.add_mesh(line1, color='red', line_width=5)
p.camera_position = [(0, 0, 1), (0.1, 0, 0), (0, 0, 0)]
# p.show()

#%%
df_lines = convert_line_df(line1, [ppu.CFD_K_SPECIES_NAME])
df_lines.plot()

#%%
dss = []
for key, fp in cfd_fp_dict.items():
    phi, kwt = key.split('_')
    mesh = pv.read(fp)
    line_out = extract_line_axi(mesh, a, b)
    df_lines = convert_line_df(line_out, cfd_all_fields)
    df_int = interp_df_to_new_index(df_lines, beam_path_dist_grid.to('m').magnitude)
    ds_out = xr.Dataset(df_int)
    ds_out = ds_out.rename({'dim_0': 'dist'})
    ds_out = ds_out.assign_coords(kwt=float(kwt))
    ds_out = ds_out.assign_coords(phi=float(phi))
    ds_out.expand_dims('kwt').expand_dims('phi').stack(temp=('kwt', 'phi'))
    dss.append(ds_out)

ds_out = xr.concat(dss, 'temp')
ds_out = ds_out.set_index(temp=['kwt', 'phi']).unstack('temp')
ds_out['dist'] = ds_out['dist'].pint.quantify('m').pint.to('mm')
ds_out.pint.dequantify().to_netcdf(pjoin('output', 'cfd_profiles_beam_barrelexit.cdf'))

#%%
ds_out[ppu.CFD_K_SPECIES_NAME].plot(col='kwt', hue='phi')

#%%[markdown]
# Extract lines along torch axis

#%%
a = [CFD_EXIT_OFFSET + Quantity(0, 'cm'), Quantity(0, 'cm'), Quantity(0, 'cm')]
b = [CFD_EXIT_OFFSET + Quantity(40, 'cm'), Quantity(0, 'cm'), Quantity(0, 'cm')]
line1 = extract_line_axi(mesh, a, b, resolution=1000)

p = pv.Plotter()
p.add_mesh(mesh, scalars='K')
p.add_mesh(line1, color='red', line_width=1)
p.camera_position = [(0, 0, 1), (0.1, 0, 0), (0, 0, 0)]
# p.show()

#%%

#TODO: add multi file for beam positions above. 

dss = []
for key, fp in cfd_fp_dict.items():
    for offset in beam_offsets:
        x_start = 0
        a = [CFD_EXIT_OFFSET + Quantity(x_start, 'cm'), Quantity(0, 'cm'), offset]
        b = [CFD_EXIT_OFFSET + Quantity(40, 'cm'), Quantity(0, 'cm'), offset]
        phi, kwt = key.split('_')
        mesh = pv.read(fp)
        line_out = extract_line_axi(mesh, a, b, resolution=1000)
        df_lines = convert_line_df(line_out, cfd_all_fields)

        # No interpolation needed here
        
        # df_int = interp_df_to_new_index(df_lines, dist_grid.to('m').magnitude)
        df_int = df_lines
        ds_out = xr.Dataset(df_int)
        ds_out = ds_out.rename({'dim_0': 'x'})
        ds_out = ds_out.assign_coords(kwt=float(kwt))
        ds_out = ds_out.assign_coords(phi=float(phi))
        ds_out = ds_out.assign_coords(offset=offset.to('m').magnitude)
        ds_out.expand_dims('kwt').expand_dims('phi').expand_dims('offset').stack(temp=('kwt', 'phi', 'offset'))
        dss.append(ds_out)

ds_out = xr.concat(dss, 'temp')
ds_out = ds_out.set_index(temp=['kwt', 'phi', 'offset']).unstack('temp')
ds_out['x'] = ds_out['x'].pint.quantify('m').pint.to('mm')
ds_out['offset'] = ds_out['offset'].pint.quantify('m').pint.to('mm')
ds_out.pint.dequantify().to_netcdf(pjoin('output', 'cfd_profiles_centerline.cdf'))

#%%
df_lines
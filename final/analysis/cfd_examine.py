#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from pi_paper_utils.constants import *

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

goldi_pos = Quantity(180, 'mm')

#%%

#%%
fp_cfd_2d = pjoin(REPO_DIR, 'final','dataset','output', 'mdot0130_phi080_K100.csv')

df_cfd = pd.read_csv(fp_cfd_2d, index_col=[0,1])

ds_cfd = df_cfd.to_xarray()

ds_cfd['pos_x'] = (ds_cfd['pos_x'] * 1000) - 208
ds_cfd['pos_y'] = ds_cfd['pos_y'] * 1000

#%%

ds_cfd['K'].plot(x = 'pos_x', robust=True)

#%%

ds_cfd['Yeq_K+'].plot(x = 'pos_x', robust=True)

#%%

ds_cfd.sel(pos_x=50, method='nearest')[['Yeq_K','K']].to_array('var').plot(hue='var')

#%%

ds_cfd.sel(pos_y=0, method='nearest')[['Yeq_K','K']].to_array('var').plot(hue='var')
plt.yscale('log')




# %%
ds_cfd_beam = ppu.fileio.load_cfd_beam(convert_rho_number=True)

ds_cfd_beam = ds_cfd_beam.sel(phi=0.8).sel(offset=0)

da_cfd_beam = ds_cfd_beam[CFD_K_SPECIES_NAME]
# da_cfd_beam = ds_cfd_beam['Yeq_K']  
da_cfd_beam = da_cfd_beam/da_cfd_beam.max('dist')

da_cfd_beam
da_cfd_beam.sel(motor=goldi_pos,method='nearest').plot(col='mp', hue='kwt')

plt.yscale('log')

# plt.ylim(1e-7,)
plt.ylim(1e-3,2)
plt.xlim(3,8)

for ax in plt.gcf().axes:
    ax.set_xlabel('Beam Distance [cm]')

plt.savefig(pjoin(DIR_FIG_OUT, 'cfd_beam_goldi_pos_kwt.png'))

#%%


ds_K = ds_cfd_beam[['Yeq_K', 'K']].to_array('var').sel(kwt=1, method='nearest')#.plot(hue='var',row='kwt')

ds_K = ds_K.sel(motor = [50, 100, 150,175,200,250], method='nearest')

ds_K.plot(hue='motor', col='var', row='mp')

plt.yscale('log')

#%%

ds_K.sel(mp='mw_horns').sel(motor=180, method='nearest').plot(hue='var')

plt.yscale('log')

plt.ylim(1e-10,)



#%%

plt.figure(figsize=(6,3))

da_cfd_pos_plot = da_cfd_beam.sel(kwt=1, method='nearest').sel(mp='mw_horns')
da_cfd_pos_plot = da_cfd_pos_plot.sel(motor = [50, 100, 150,175,200,250], method='nearest')
da_cfd_pos_plot = da_cfd_pos_plot.assign_coords(motor = da_cfd_pos_plot.coords['motor'].round(2))

da_cfd_pos_plot.plot(hue='motor')

plt.yscale('log')

plt.ylim(1e-7,)
plt.xlim(3,8)

plt.gca().get_legend().set_bbox_to_anchor((1.05,0.5))

plt.gca().set_xlabel('Beam Distance [cm]')


plt.savefig(pjoin(DIR_FIG_OUT, 'cfd_beam_pos_dep.png'))

# %%
ds_cfd = ppu.fileio.load_cfd_centerline()

ds_cfd = ds_cfd.sel(offset=0).sel(phi=0.8)

ds_cfd['nK_m3'] = ds_cfd['Yeq_K'].pint.to('particle/m^3')

#%%


all_K_species = ['Yeq_K','Yeq_K+','Yeq_K2CO3','Yeq_KO','Yeq_KOH']
ds_cfd[[*all_K_species, 'all_K_Yeq']].to_array('var').plot(hue='var',row='kwt')

plt.yscale('log')

# plt.gca().get_legend().set_bbox_to_anchor((1,1))

plt.ylim(1e8,1e17)

# plt.savefig(pjoin(DIR_FIG_OUT, 'cfd_K_species.png'), dpi=300, bbox_inches='tight')



# %%

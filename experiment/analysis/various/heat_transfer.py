#%%
from mhdlab.analysis.standard_import import *
from mhdlab.fileio.ct import load_df_cuttimes
from mhdlab.coords import assign_signal, unstack_multindexed_acq_dim

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

fp_dsst = pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')
dsst = mhdlab.fileio.TFxr(fp_dsst).as_dsst(convert_to_PT=False)

fp_dst_coords = pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdlab.fileio.TFxr(fp_dst_coords).as_dsst(convert_to_PT=False)['coords']

dst_coords

#%%


df = pd.read_csv(pjoin(REPO_DIR, 'modeling','cfd','simulation_data', 'heatFlux_last.csv'))

s = df.set_index(['phi', 'Y_K'])['inner_total']

da_cfd = xr.DataArray(s).unstack('dim_0')
da_cfd = da_cfd*(-1)
da_cfd = da_cfd.rename(Y_K='kwt')
da_cfd.coords['kwt'] = da_cfd.coords['kwt']*100

da_cfd.plot(hue='kwt')

# %%
# %%[markdown]

# # Seed Ramp

#%%

fp_cuttimes = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_testcase_kwt.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes).sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')

df_cuttimes['date'] = df_cuttimes['Start Time'].dt.date

#%%

hts = []

for idx, row in df_cuttimes.iterrows():
    timeslice = slice(row['Start Time'], row['Stop Time'])

    ht = dsst['calorimetry']['CC_heatTransfer'].sel(time=timeslice)
    ht = assign_signal(ht, dst_coords['kwt'].round(2), timeindex='time')

    hts.append(ht)

ht = xr.concat(hts, dim='time')

ht = unstack_multindexed_acq_dim(ht)

#%%

ht['CC_heatTransfer'].mean('mnum').plot(marker='o', label='expt')

da_cfd.sel(phi=0.8).plot(label='cfd', marker='o')

plt.legend()

#%%

from mhdlab.plot import xr_errorbar

mean = ht['CC_heatTransfer'].mean('mnum')
std = ht['CC_heatTransfer'].std('mnum')

xr_errorbar(mean, std)
# %%[markdown]

# # Equivalence Ratio


#%%
fp_cuttimes = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_testcase_phi.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes).sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')

df_cuttimes['date'] = df_cuttimes['Start Time'].dt.date
# %%
hts = []

for idx, row in df_cuttimes.iterrows():
    timeslice = slice(row['Start Time'], row['Stop Time'])

    ht = dsst['calorimetry']['CC_heatTransfer'].sel(time=timeslice)
    ht = assign_signal(ht, dst_coords['phi'].round(2), timeindex='time')
    ht = assign_signal(ht, dst_coords['kwt'].round(1), timeindex='time')

    hts.append(ht)

ht = xr.concat(hts, dim='time')

ht = unstack_multindexed_acq_dim(ht)

ht = ht['CC_heatTransfer']
#%%

ht.mean('mnum').plot(hue='kwt', marker='o')

#%%

mean = ht.mean('mnum')
std = ht.std('mnum')

xr_errorbar(mean, std, huedim='kwt')



#%%

fig, axes = plt.subplots(1, 2, figsize=(10,5), sharey=True)

ht.mean('mnum').plot(hue='kwt', marker='o', ax=axes[0])
axes[0].set_title('Experiment')
da_cfd.plot(hue='kwt', ax=axes[1])
axes[1].set_title('CFD')



# %%

ds = xr.merge([
    ht.mean('mnum').rename('expt'),
    da_cfd.rename('cfd')
])


fig, ax = plt.subplots(1, 1, figsize=(5,5))


colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

lns_expt = []
lns_cfd = []

for kwt in ds['kwt']:
    color = colors.pop(0)
    ln_expt = ds['expt'].dropna('phi', 'all').sel(kwt=kwt).plot(x='phi', ax=ax, marker='o', label=f'kwt={kwt}', color=color, linestyle='-')
    ln_cfd = ds['cfd'].dropna('phi','all').sel(kwt=kwt).plot(x='phi', ax=ax, marker='o', label=f'kwt={kwt}', linestyle='--', color=color)

    lns_expt.append(ln_expt[0])
    lns_cfd.append(ln_cfd[0])

lns = lns_expt + lns_cfd

labs = ['expt kwt=0%', 'expt kwt=0.1%', 'expt kwt=1%', 'cfd kwt=0%', 'cfd kwt=0.1', 'cfd kwt=1%']
ax.legend(lns, labs)

plt.ylabel('Heat Transfer (kW)')
plt.xlabel('Equivalence Ratio')

plt.savefig(pjoin(DIR_FIG_OUT, 'heat_transfer_phi_cfd_expt.png'))


#%%

lns
    
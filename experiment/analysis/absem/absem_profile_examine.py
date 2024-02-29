#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis.absem.fitting import Q_2peak




# %%

Ts = np.logspace(0,4,5)

wls = np.arange(760,775,0.001)


for T in Ts:
    T = Quantity(T, 'K')
    Qs = Q_2peak(wls, T=T)
    plt.plot(wls, Qs, label=T)


plt.legend()
# plt.yscale('log')

plt.xlim(766.4,766.6)
# %%

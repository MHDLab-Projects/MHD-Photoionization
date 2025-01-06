#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis.absem.fitting import kappa_2peak




# %%

Ts = np.logspace(0,4,5)

wls = np.arange(760,775,0.001)


for T in Ts:
    T = Quantity(T, 'K')
    kappas = kappa_2peak(wls, T=T)
    plt.plot(wls, kappas, label=T)


plt.legend()
# plt.yscale('log')

plt.xlim(766.4,766.6)
# %%

# Calculating collisional broadening from Alkemade II.231 and Table VII.9


# sigma = Quantity(60.4, 'angstrom^2')
sigma = Quantity(57.7, 'angstrom^2')
T = Quantity(2500, 'K')

m_N2 = Quantity(28, 'amu')
m_K = Quantity(39, 'amu')

redmass_N2 = m_N2*m_K/(m_N2+m_K)

kb = Quantity(1, 'boltzmann_constant')

p = Quantity(1, 'atm')

nx = p/(kb*T)

vth = np.sqrt(8*kb*T/(redmass_N2*np.pi))

# Alkemade formula. Not sure where this pi is coming from
# nu_c = vth*sigma*nx/np.pi

nu_c = vth*sigma*nx

nu_c.to('1/ns')
#%%

vth.to('cm/s')

#%%

nx.to('1/m^3')
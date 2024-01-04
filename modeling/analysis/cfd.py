#%%

from mhdpy.analysis.standard_import import *

import pint_pandas
from pint import Quantity

fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles.csv')

df = pd.read_csv(fp, index_col=0)

df['T'] = df['T'].astype('pint[K]')
df['p'] = df['p'].astype('pint[Pa]')

#%%

R = Quantity(8.31446261815324, 'J/(mol*K)')
df['rho'] = df['p']/(df['T']*R)

NA = Quantity(6.02214076e23, '1/mol')

df['rho'] = df['rho']*NA

df['rho'].pint.to('1/cm^3')
#%%

df['rho'].pint.magnitude.plot()


#%%

Kp_rho = df['Kp']*df['rho']
Kp_rho = Kp_rho.pint.to('1/cm^3')

Kp_rho
#%%

Kp_rho.pint.magnitude.plot()

plt.yscale('log')
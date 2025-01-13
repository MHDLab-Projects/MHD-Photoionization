"""Collect recipes from test plan and save them to a file. """
#%%

import os
from mhdlab.analysis.standard_import import *

munged_dir = pjoin(REPO_DIR, 'experiment', 'data', 'munged')

ss = []

for date in ['2023-04-07', '2023-05-12','2023-05-18', '2023-05-24']:
    fp = pjoin(munged_dir, date, 'TestPlan.xlsx')
    df = pd.read_excel(fp, sheet_name='Emulsion Recipe', index_col=0)
    s = df['values']
    s.name = date
    ss.append(s)

df = pd.concat(ss, axis=1)


#%%

# As noted in 2023-05-24, the fuel amount was slightly off in the other test plans.
# This was corrected in the 2023-05-24 test plan from 109 g to 106.5g 

recipe_sel = df[['2023-05-18', '2023-05-24']]

recipe_sel.rename({'2023-05-18': '2023-04-07 to 2023-05-18'}, axis=1, inplace=True)

recipe_sel.to_csv(pjoin(DIR_DATA_OUT, 'em_recipe_2023-05-24.csv'))


# %%


"""
Generate presentation table from sim_input_mean.xlsx 
"""

#%%

import pandas as pd
import numpy as np
import os
from dotenv import load_dotenv
load_dotenv()

from mhdlab.analysis.standard_import import *


fp = pjoin('output', 'sim_input_mean.xlsx')

sheet_names = ['HVOF Process Inputs', 'Calorimetry', 'Other Temperatures']

dfs = pd.read_excel(fp, sheet_name=sheet_names, header=[0,1,2], index_col=[0,1,2])

#%%

df_hvof = dfs['HVOF Process Inputs']
df_cal = dfs['Calorimetry']

df = pd.concat([df_hvof, df_cal], axis=1)

df = df.dropna(how='all')   

# Hacky wat to convert mass fraction to percent
df['CC_K_massFrac_in'] = df['CC_K_massFrac_in']*100

# Get the current column MultiIndex
columns = df.columns

# Create a list to hold the new column names
new_columns = []

# Loop over the current column names
for col in columns:
    if col[0] == 'CC_K_massFrac_in':
        # If the first level of the column name is 'CC_K_massPct_in', change the second and third levels
        new_columns.append(('CC_K_massPct_in', '%', 'K Mass Percent'))
    else:
        # Otherwise, keep the column name as it is
        new_columns.append(col)

# Replace the current column MultiIndex with the new column names
df.columns = pd.MultiIndex.from_tuples(new_columns)




df.info()



#%%

downselct_cols = [
    'CC_total_flow_in', 
    'CC_equivalenceRatio',
    'CC_K_massPct_in',
    'CC_heatInput',
    'CC_heatTransfer'
    ]


# Downselect columns
df = df.loc[:, df.columns.get_level_values(0).isin(downselct_cols)]

df

#%%


#%%

# Split the DataFrame into two based on the first level of the index
df_mean = df.xs('mean', level=0)
df_std = df.xs('std', level=0)



#%%

df_mean
# %%
# Create a new DataFrame where each cell is a string that contains the mean and standard deviation values
df_combined = df_mean.applymap('{:.2f}'.format) + " ± " + df_std.applymap('{:.2f}'.format)

# Get the current column names
columns = df_combined.columns

# Create new column names by combining the second and third levels
new_columns = [f"{long_name} [{unit}]" for long_name, unit in zip(columns.get_level_values(2), columns.get_level_values(1))]

# Replace the current column MultiIndex with the new column names
df_combined.columns = new_columns
# %%


df_combined
# %%

df_combined.to_csv('output/sim_input_mean_table.csv')

#%%

df_combined.to_clipboard()
# %%

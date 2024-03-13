#%%

# Generates a report for the simulation input

#%%

from mhdpy.analysis.standard_import import *


# load and excel file

# load the excel file
import pandas as pd


# Replace 'filename.xlsx' with the path to your Excel file
df_hvof = pd.read_excel(pjoin('output', 'sim_input_all.xlsx'), sheet_name='HVOF Process Inputs', header=[0,1,2], index_col=[0,1])
df_calor = pd.read_excel(pjoin('output', 'sim_input_all.xlsx'), sheet_name='Calorimetry', header=[0,1,2], index_col=[0,1])

df = pd.concat([df_hvof, df_calor], axis=1)

units = pd.Series(df.columns.get_level_values(1), index=df.columns.get_level_values(0))
long_names = pd.Series(df.columns.get_level_values(2), index=df.columns.get_level_values(0))

# Drop unit and long name information
df.columns = df.columns.droplevel(2)
df.columns = df.columns.droplevel(1)


df = df[['CC_em_flow_in', 'CC_fuel_flow_in', 'CC_o2_flow_in', 'CC_P', 'CC_water_flow_in', 'CC_water_T_in', 'CC_water_T_out']]

#%%


df.columns

#%%


df = df.reset_index()

 # Reset the index level you want to split

# Split the reset index
df[['kwt', 'date', 'num']] = df['Test Case'].str.split('_', expand=True)

df['repeat'] = df['date'] + '_' + df['num']

df = df.sort_values(['repeat'])

# Set the new MultiIndex
df.set_index(['kwt', 'repeat', 'Statistic'], inplace=True)

df = df.drop(['Test Case', 'date', 'num'], axis=1)

#%%
from pylatex import Document, Section, Figure, NoEscape

output_dir = pjoin(DIR_FIG_OUT, 'signal_stats_repeat')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for col in df.columns:
    df_col = df[col]


    kwt_gb = df_col.groupby('kwt')

    for kwt, df_kwt in kwt_gb:
        fig, ax = plt.subplots()

        repeat_values = df_col.index.get_level_values('repeat').unique()
        plt.xticks(ticks=range(len(repeat_values)), labels=repeat_values, rotation=90)
        plt.xlabel('Repeat')

        df_kwt = df_kwt.droplevel('kwt')

        s = df_kwt

        s_mean = s.xs('mean', level='Statistic')
        s_std = s.xs('std', level='Statistic')

        # Plot with error bars
        plt.errorbar(x=s_mean.index, y=s_mean, yerr=s_std, capsize=4, marker='o', label=kwt)

        # Set x-ticks to all unique values in the index

        plt.ylabel("{} ({})".format(col, units[col]))
        plt.title(long_names[col])


        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='kwt')
        plt.tight_layout()


        plt.savefig(pjoin(output_dir, '{}_{}.png'.format(kwt.replace('.', 'p'), col)))
        plt.close()

# doc.generate_pdf(clean_tex=False)

# %%

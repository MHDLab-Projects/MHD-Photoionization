#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
from mhdpy.analysis import mws
import seaborn as sns
import pi_paper_utils as ppu
# sns.set_theme(style="darkgrid")

def main(datestr):

    tcs = tc_dict[datestr]

    data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
    DIR_LECROY_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

    dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

    figure_out_dir = pjoin(DIR_DATA_OUT, 'mws_1d', datestr)

    dss = {}
    dss_std = {}

    for tc in tcs:

        fp_in = pjoin(DIR_LECROY_PROC_DATA, '{}.cdf'.format(tc))

        mag_0 = ppu.fileio.load_mws_T0()

        ds_in = xr.load_dataset(fp_in)

        ds_in = ds_in.mws.calc_AS_abs(mag_0=mag_0).mws.calc_time_stats()
        ds_in = ds_in.sel(date=datestr).sel(run_num=1)
        ds_in = ds_in[['i','q','mag','phase','AS_abs', 'dAS_abs']] #Plotted variables

        tc_dim = [dim for dim in ds_in.dims if dim not in ['time','mnum']][0]
        mnum_counts = ds_in['i'].mean('time').groupby(tc_dim).count('mnum')

        ds = ds_in.mean('mnum', keep_attrs=True)
        ds_std = ds_in.std('mnum', keep_attrs=True)

        tc_fig_dir = pjoin(figure_out_dir, tc)
        if not os.path.exists(tc_fig_dir): os.makedirs(tc_fig_dir)

        for fn in os.listdir(tc_fig_dir):
            os.remove(os.path.join(tc_fig_dir, fn))

        plt.figure()

        mnum_counts.plot(marker='o')
        plt.ylabel('Acquisitions')


        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'acquisition_counts.png'))


        dss[tc] = ds
        dss_std[tc] = ds_std

    # %%

    figsize = (8, 11)
    for tc in dss:
        ds = dss[tc]

        tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]

        plot_names = np.array([ds.data_vars]).T
        fig, axes = plt.subplot_mosaic(plot_names, figsize=figsize, sharex=True)

        for var in ds.data_vars:
            ax = axes[var]
            ds[var].plot(ax=ax, hue=tc_dim)

            if var != 'i': ax.get_legend().remove()
            if var != 'i': ax.set_title('')

        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'raw_signals_hue.png'))

    # %%

    figsize = (8, 11)
    for tc in dss:
        ds = dss[tc]

        plot_names = np.array([ds.data_vars]).T
        fig, axes = plt.subplot_mosaic(plot_names, figsize=figsize, sharex=True)

        for var in ds.data_vars:
            ax = axes[var]
            ds[var].plot(ax=ax)

        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'raw_signals_color.png'))


    #%%

    for tc in dss:

        plt.figure()
        ds = dss[tc]

        tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]

        ds['AS_abs'].plot(hue=tc_dim)

        plt.yscale('log')

        plt.xlim(-1,30)

        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'AS_log.png'))
    #%%


    for tc in dss:

        ds = dss[tc]
        ds_std = dss_std[tc]

        tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]

        plt.figure()

        tc_vals = ds.coords[tc_dim]
        fig, axes = plt.subplots(len(tc_vals), figsize=(5, 3*len(tc_vals)), sharex=True)
        # ds['AS_abs'].plot(hue=tc_dim)
        # dss_std[tc]['AS_abs'].plot(hue=tc_dim, linestyle='--')

        ds = ds.drop_vars([var for var in ds.coords if var not in ds.dims])

        # xr_errorbar(ds['AS_abs'], dss_std[tc]['AS_abs'], huedim=tc_dim)

        for i, c in enumerate(ds.coords[tc_dim]):
            ds_sel_mean = ds.sel({tc_dim: c})['AS_abs']
            xs = ds_sel_mean.coords['time'].values

            vals_mean = ds_sel_mean.values
            vals_std = ds_std.sel({tc_dim: c})['AS_abs'].values

            plt.sca(axes[i])
            plt.plot(xs, vals_mean)

            plt.fill_between(x=xs,
                    y1=vals_mean - vals_std,
                    y2=vals_mean + vals_std,
                    alpha=0.5
                    )


        # plt.yscale('log')

        # plt.xlim(-1,30)
        # plt.ylim(-1e-2, 1e-2)


        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'AS_lin.png'))

    #%%


    # Load an example dataset with long-form data
    fmri = sns.load_dataset("fmri")

    # Plot the responses for different events and regions
    sns.lineplot(x="timepoint", y="signal",
                hue="region", style="event",
                data=fmri)


# %%

    # %%

    max_window = slice(-1,1)
    pre_pulse = slice(-50,-1)

    for tc in dss:
        ds = dss[tc]

        plot_names = [['pulse_max'], ['pre_pulse_avg'], ['pre_pulse_std']]
        fig, axes = plt.subplot_mosaic(plot_names, figsize=figsize, sharex=True)

        pulse_max = ds['AS_abs'].sel(time=max_window).max('time')
        pulse_max.name = 'pulse_max'
        pre_pulse_avg =ds['AS_abs'].sel(time=pre_pulse).mean('time') 
        pre_pulse_avg.name = 'pre_pulse_avg'
        pre_pulse_std =ds['AS_abs'].sel(time=pre_pulse).std('time') 
        pre_pulse_std.name = 'pre_pulse_std'

        ds_out = xr.merge([
            pulse_max,
            pre_pulse_avg,
            pre_pulse_std
        ])

        for var in ds_out.data_vars:
            ax = axes[var]
            ds_out[var].plot(ax=ax, marker='o')
        
        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'pulse_stats.png'))


    # %%

    from mhdpy.analysis.mws.fitting import pipe_fit_mws_1 

    das_fit_input = {}
    dss_fits = {}
    dss_p = {}
    dss_p_stderr ={}

    for tc in dss:
        print(tc)
        ds = dss[tc]

        da_fit = ds['dAS_abs'].copy()

        ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_1(da_fit)

        das_fit_input[tc] = da_fit
        dss_fits[tc] = ds_mws_fit['AS_fit']
        dss_p[tc] = ds_p
        dss_p_stderr[tc] = ds_p_stderr

        # break
    # %%
    for tc in dss_p:

        plt.figure()
        ds_p = dss_p[tc]

        tc_dim = [dim for dim in ds_p.dims if dim != 'time'][0]

        g = ds_p.to_array('var').plot(row='var', yscale = 'log', sharey=False, marker='o')

        # g.axes[1].set_ylim()


        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'ds_p.png'))

    #%%

    for tc in dss_fits:
        plt.figure()

        fits = dss_fits[tc]
        da_fit = das_fit_input[tc]


        tc_dim = [dim for dim in fits.dims if dim != 'time'][0]
        tc_vals = da_fit.coords[tc_dim]

        fig, axes = plt.subplots(len(tc_vals), figsize=(5, 3*len(tc_vals)), sharex=True)

        for i, tc_val in enumerate(tc_vals):

            ax = axes[i]

            da_fit.sel({tc_dim: tc_val}).plot(ax=ax)
            fits.sel({tc_dim: tc_val}).plot(ax=ax, linestyle='--', color='gray')

            ax.set_yscale('log')
            # plt.ylim(-0.1,1.1)
            ax.set_xlim(-1,30)

            if i != len(tc_vals) -1: ax.set_xlabel('')

            # if i != 0:
            #     ax.set_title('')

        tc_fig_dir = pjoin(figure_out_dir, tc)
        plt.savefig(pjoin(tc_fig_dir, 'fits_all.png'))


tc_dict = {
    '2023-04-07': [
        '536_pos',
        '536_power',
        '53x',
    ],
    '2023-05-12': [
        '536_pos',
        '536_power',
        '53x',
    ],
    '2023-05-18': [
        '536_pos',
        '536_power',
        '53x',
    ],
    '2023-05-24': [
        # '516_pos_1',
        '536_power',
        '53x'
    ],
}

if __name__ == '__main__':
    for datestr in tc_dict:
        print(datestr)
        main(datestr)